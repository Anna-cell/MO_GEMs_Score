#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 14:57:56 2023

@author: e158401a
"""

import pandas as pd
import cobra
import mocbapy
from mocbapy.EcosystemModel import create_model, bensolve_default_options
import mocbapy.analysis
import diet
import pickle
import numpy as np

diet_param = "Western_diet"
host_choice = "sIEC"
semiconstrained = False
block = False
sampling = False
save = False

link_GCA_taxo = pd.read_csv("resource/link_taxo_GCA_ncbi_16bact.csv", index_col="GCA")

host = cobra.io.read_sbml_model("resources/enterocyte_ori.xml")
host.objective = host.reactions.get_by_id('biomass_reactionIEC01b')
host.solver = 'cplex'

all_metabolites = []
for met in host.metabolites:
    all_metabolites.append((met.id, host))
host_int = diet.medium_2(model = host.copy(), medium_lumen = diet_param, host_choice = host_choice,
                         block = block, mocba=False, semiconstrained = semiconstrained)[0]
growth_host = host_int.optimize().objective_value
host, diet_dict = diet.medium_2(model = host, medium_lumen = diet_param, host_choice = host_choice, 
                                  block = block, mocba=True, semiconstrained = semiconstrained)

#bact_selection = [0, 1, 3, 4, 5, 6, 8, 12, 13, 15]#selection of ten bacteria with less niche overlap (see Figure S5)
bact_selection = [0, 2, 9, 11] #In order : A. muciniphila, B. xylanislvens, F. prausnitzii, R. bromii


#Normal ecosystem pareto
bacteria_models = []
growth_b = []
bacteria_ids = []
for bact in link_GCA_taxo.index[bact_selection]:
    b = cobra.io.read_sbml_model("resources/minimal_microbiome_models/"+bact+"_carve.xml")
    b.solver = "cplex"
    bacteria_id = link_GCA_taxo.loc[bact,"taxonomy"]
    print(bacteria_id)
    bacteria_ids.append(bacteria_id.replace(" ","_"))
    bacterimet = []
    for met in b.metabolites:
        bacterimet.append((met.id, b))
    all_metabolites = list(set(all_metabolites + bacterimet))
    b_int = diet.medium_2(model = b.copy(), medium_lumen = diet_param, host_choice = "bacteria", 
                          block = block, mocba=False, semiconstrained = semiconstrained)[0]
    growth_b.append(b_int.optimize().objective_value)
    b, diet_dict_b = diet.medium_2(model = b, medium_lumen = diet_param, host_choice = "bacteria", 
                                   block=block, mocba=True, semiconstrained = semiconstrained)
    diet_dict = {**diet_dict_b, **diet_dict} #When a key is the same, only the value from the second dico stays. 
    bacteria_models.append(b)

diet_dict["o2_e"] = (0, host.reactions.get_by_id("EX_o2_e").upper_bound)

metabolic_dict = {x:x[0] for x in all_metabolites}
ecosys = create_model(model_array=[host]+bacteria_models, metabolic_dict=metabolic_dict, diet = diet_dict)

bensolve_opts = bensolve_default_options()
bensolve_opts['message_level'] = 0
sol_mofba = mocbapy.analysis.mo_fba(ecosys, options=bensolve_opts)

ids = ["enterocyte"]+bacteria_ids
growth_ini = [growth_host]+growth_b
points = sol_mofba.Primal.vertex_value

## To save
if save:
    ids = ["enterocyte"]+bacteria_ids
    growth_ini = [growth_host]+growth_b
    points = sol_mofba.Primal.vertex_value
    with open("sol_mofba_minimal_microbiome_Pareto_extremes.pickle", "wb") as fp:
        pickle.dump([sol_mofba.Primal.vertex_type, sol_mofba.Primal.vertex_value, ids, growth_ini], fp)

# To analyse
growth_alone_eco = pd.DataFrame({"ids" : ids, "growth_ini" : growth_ini})
for i in range(len(sol_mofba.Primal.vertex_type)):
    if  sol_mofba.Primal.vertex_type[i]:
        growth_alone_eco[i] =  sol_mofba.Primal.vertex_value[i]

growth_alone_eco.index = growth_alone_eco["ids"]
growth_alone_eco.drop("ids", inplace=True, axis = 1)
growth_alone_eco_T = growth_alone_eco.transpose()
growth_alone_eco_T["Ecosystem_biomass"] = list(growth_alone_eco.sum())

# sampling (for suboptimal points on the PCA), need computational power !
if sampling:
    # First, translate mocbapy object into cobra model
    model = cobra.Model('ecosys')
    for m in range(len(ecosys.sysmetabolites)):
        model.add_metabolites(cobra.Metabolite(ecosys.sysmetabolites[m]))
    for r in range(len(ecosys.sysreactions)):
        reaction = cobra.Reaction(ecosys.sysreactions[r])
        reaction.lower_bound = ecosys.lb[r]
        reaction.upper_bound = ecosys.ub[r]
        dict_metabolites = {}
        model.add_reactions([reaction])
        for m in range(len(ecosys.sysmetabolites)):
            if ecosys.Ssigma[m,r] != 0:
                dict_metabolites[ecosys.sysmetabolites[m]] = ecosys.Ssigma[m,r]
        reaction.add_metabolites(dict_metabolites)
    model.solver = "cplex"
    
    #optgsampler options
    spl_size = 1000
    thinning = 100 #init = 100
    processes = None #init = None
    s = cobra.sampling.sample(model, spl_size, thinning = thinning, processes=processes)
    
    with open("sampling_ecosys.pickle", "wb") as fp:
        pickle.dump(s, fp)
 
### ANALYSIS ###   
#Calculate a score based on the hypervolume of the pareto linked to the origin point

#get rid of slutions with only -1 for one organism and 0 for all the others
points = points[np.all(points >= 0, axis=1)] 
#Add a point for the original growth of each organism, where every other organism is at 0 (maximal growth when alone)
for i in range(len(growth_ini)):
    ori = [0]*len(growth_ini)
    ori[i] = growth_ini[i]
    points = np.append(points, [ori], axis=0)
#And the original point (every growth to 0) in order to have a comparable form for each ecosystem. 
ori = [0]*len(growth_ini)
points = np.append(points, [ori], axis=0)

# Normalize each dimension by its organism's initial growth (maximal growth when alone)
for i in range(len(growth_ini)):
    points[:,i] = points[:,i] / growth_ini[i]

#Measure the hypervolume based on the given polyhedron
from scipy.spatial import ConvexHull
volume = ConvexHull(points=points)#, qhull_options='FS')
print("The ecosystem interaction score is ",volume.volume)

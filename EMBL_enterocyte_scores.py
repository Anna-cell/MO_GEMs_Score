#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 11:11:15 2022

@author: e158401a
"""
import pandas as pd
import cobra
import matplotlib.pyplot as plt
from sklearn import metrics
import mocbapy
from mocbapy.EcosystemModel import create_model, bensolve_default_options
import mocbapy.analysis
import pickle
import numpy as np
import concurrent.futures 
import diet


def pareto(b, diet_param, host, bacteria_id):
    bacterimet = []
    for met in b.metabolites:
        bacterimet.append((met.id, b))
    all_metabolites = list(set(enteromet + bacterimet))
    metabolic_dict = {x:x[0] for x in all_metabolites}
    if diet_param in ["Western_diet", "Protein_diet"]:
        host_int = diet.medium_2(model = host.copy(), medium_lumen = diet_param, host_choice = host_choice, block = False, mocba=False, semiconstrained = semiconstrained)[0]
        growth_host = host_int.optimize().objective_value
        host, diet_dict_1 = diet.medium_2(model = host, medium_lumen = diet_param, host_choice = host_choice, block = False, mocba=True, semiconstrained = semiconstrained)
        b_int = diet.medium_2(model = b.copy(), medium_lumen = diet_param, host_choice = "bacteria", block = False, mocba=False, semiconstrained = semiconstrained)[0]
        growth_b = b_int.optimize().objective_value
        b, diet_dict_2 = diet.medium_2(model = b, medium_lumen = diet_param, host_choice = "bacteria", block=False, mocba=True, semiconstrained = semiconstrained)
        diet_dict = {**diet_dict_1, **diet_dict_2}
        if o2 == "present":        
            diet_dict["o2_e"] = (-10, 1000)#host.reactions.get_by_id("EX_o2_e").upper_bound)
            #print(diet_dict["o2_e"])
        elif o2 == "absent":
            diet_dict["o2_e"] = (0, 1000)#host.reactions.get_by_id("EX_o2_e").upper_bound)
            #host.reactions.get_by_id("EX_o2_e").bounds = (-10,10) #Constrain what the enterocyte gives to the bacteria
    elif diet_param == "open":
        growth_host = host.optimize().objective_value
        growth_b = b.optimize().objective_value
        diet_dict = {}
        for ex in host.exchanges:
            met_e, suffixe = diet.no_compartment_id(list(ex.metabolites.keys())[0].id)
            diet_dict[met_e+suffixe] = ex.bounds #We keep the original bounds, that are not necessarily -1000, 1000
            ex.bounds = (-1000, 1000) #But make this "internal exchange reaction" free
        for ex in b.exchanges:
            met_e, suffixe = diet.no_compartment_id(list(ex.metabolites.keys())[0].id)
            if met_e+suffixe not in list(diet_dict.keys()): 
                diet_dict[met_e+suffixe] = ex.bounds
            else:
                bounds = diet_dict[met_e+suffixe]
                diet_dict[met_e+suffixe] = (min(bounds[0], ex.bounds[0]), max(bounds[1], ex.bounds[1])) #Keep the least restrictive bound ? Keep that in mind it might change
            ex.bounds = (-1000, 1000) #Leave the internal exchange reaction free
    ecosys = create_model(model_array=[host, b], metabolic_dict=metabolic_dict, diet = diet_dict)
    #make the fba output readable
    bensolve_opts = bensolve_default_options()
    bensolve_opts['message_level'] = 0
    sol_mofba = mocbapy.analysis.mo_fba(ecosys, options=bensolve_opts)
    points=[]
    x = []
    y = []
    #translate fba output in pareto points
    maxi_host = (growth_host, 0)
    maxi_b = (0, growth_b)
    category = ["0","0","0"]
    for i in range(len(sol_mofba.Primal.vertex_value)):
        if sol_mofba.Primal.vertex_type[i] == 1: #1 means point when 0 is a vector
            point = (sol_mofba.Primal.vertex_value[i][0], sol_mofba.Primal.vertex_value[i][1])
            points.append(point)
            x.append(sol_mofba.Primal.vertex_value[i][0]/growth_host)
            y.append(sol_mofba.Primal.vertex_value[i][1]/growth_b) 
            if sol_mofba.Primal.vertex_value[i][0] > maxi_host[0]: #Is the growth of the enterocyte at this point on the pareto higher than the previously saved maximal growth ?
                maxi_host = sol_mofba.Primal.vertex_value[i] #Then it's our new max
            if sol_mofba.Primal.vertex_value[i][1] > maxi_b[1]:
                maxi_b = sol_mofba.Primal.vertex_value[i]
    if maxi_host[0] > growth_host+(0.01 * growth_host): #growth + 1% of itself to avoid approximations. 
        added_growth_host = maxi_host[0]-growth_host
        category[0] = "1"
    else:
        added_growth_host = 0
    if maxi_b[1] > growth_b+(0.01 * growth_b): #growth + 1% of itself to avoid approximations. 
        added_growth_b = maxi_b[1]-growth_b
        category[1] = "1"
    else:
        added_growth_b = 0
    if tuple(maxi_host) == tuple(maxi_b) and (category[0] == "1" and category[1] == "1"):
        category[2] = "1"
    category=''.join(category)
    nb_points = len(x)
    #Add extreme points to make a full pareto, not only the variable front
    xy = pd.DataFrame({'x': x, 'y': y})
    xy.sort_values('x', inplace=True)
    if not ((xy['x'] == 1) & (xy['y'] == 0)).any():
        xy = xy.append({'x' : 1.00001, 'y' : -0.00001}, ignore_index=True) #Both are slightky dfferent than 1 and 0 because when a serie of coordinate is not continuous it's a problem
        #xy = xy.append({'x': (growth_host/norm_host)+(norm_host/10000000000), "y":0}, ignore_index=True) #because if the points are equal it's bad. Thi one is for there is no norm, so to work on maybe
    if not ((xy['x'] == 0) & (xy['y'] == 1)).any():
        xy.loc[-1] = [-0.00001, 1.00001]
        xy.index = xy.index + 1  # shifting index
        xy = xy.sort_index()  # sorting by index
    try:
        AUC = metrics.auc(x = xy['x'], y = xy['y'])
    except: #If the growth alone is lower than with bacteria, making it hard to calculate AUC (x is not monotonous since it goes up and then down)
    #In this case, I calculate the area under the curve where x is monotonous, and the last part of the curve come back and I have to take back the value of it's AUC.    
        AUC1 = metrics.auc(x = xy['x'][0:-1], y = xy['y'][0:-1])
        AUC2 = metrics.auc(x = xy['x'][-2:], y = xy['y'][-2:])
        AUC = AUC1 - AUC2
    if AUC < 0.999 and category == "000":
        category = "-000"
    elif AUC > 0.999 and AUC < 1.0001:
        category = "=000"    
    if plot:
        plt.title("AUC = "+str(AUC)+"\n ori_b ="+str(growth_b) +" bwe = "+str(max(y)*growth_b)+"\n category =" + category)
        plt.xlabel(host.id+"'s growth")
        plt.ylabel(bacteria_id+"'s growth")
        plt.plot(xy['x'], xy['y'], 'r', linestyle="-")
        plt.axhline(y = 1, color = 'b', linestyle = '--')
        plt.axvline(x = 1, color = 'b', linestyle = '--')
        plt.show()
    return AUC, nb_points, added_growth_host, added_growth_b, growth_b, category, len(b.exchanges)

diet_param = "Western_diet"
semiconstrained = False
host_choice = "sIEC"
plot = True
o2 = "absent"

emblpath = "resources/embl_gems/"
epath = "resources/enterocyte_ori.xml" 

host = cobra.io.read_sbml_model(epath)
host.objective = host.reactions.get_by_id('biomass_reactionIEC01b')
host.solver = 'cplex'

enteromet = []
for met in host.metabolites:
    enteromet.append((met.id, host))

carve_models = pd.read_csv(emblpath+"/model_list.tsv", sep = "\t", header = 0, index_col=1)
with open("resources/carve_in_vmh.pickle", "rb") as fp:   # Unpickling
    carve_in_vmh = pickle.load(fp)

def apply_pareto_carveinVMH(i):
    model = np.int64(carve_in_vmh[i]) # If i'm iterating
    #model = np.int64(i) # If I have ncbi id
    path = carve_models.loc[model, "file_path"][0:-3]
    b = cobra.io.read_sbml_model(emblpath+path)
    bacteria_id = carve_models.loc[model, "organism_name"]
    AUC, nb_point, added_growth_host, added_growth_b, growth_b, category, nb_b_exchanges = pareto(b= b, diet_param = diet_param, host = host, bacteria_id = bacteria_id)
    return bacteria_id, AUC, nb_point, added_growth_host, added_growth_b, i, growth_b, category, nb_b_exchanges

pool_bacteria = carve_in_vmh #Must be a list of model accession in this case. But the pareto can be called with another code the apply_pareto and from the model + bacteria_id. 

if __name__ == '__main__':
    results = pd.DataFrame(columns=["ncbi_id", "bacteria_id", "AUC", "nb_points","added_growth_host", "added_growth_b", "growth_b", "category", "nb_b_exchanges"])
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        for result in executor.map(apply_pareto_carveinVMH, range(len(pool_bacteria))): #If iterating
            results = results.append({"ncbi_id":carve_in_vmh[result[5]], "bacteria_id":result[0] ,"AUC":result[1], "nb_points":result[2], "added_growth_host":result[3], "added_growth_b":result[4], "growth_b":result[6], "category" : result[7], "nb_b_exchanges" : result[8]}, ignore_index=True)

results.set_index(results["ncbi_id"])

#with open("full_info_df_"+diet_param+"_"+host_choice+"_semiconstrained"+str(semiconstrained)+".pickle", "wb") as fp:
#    pickle.dump(results, fp)

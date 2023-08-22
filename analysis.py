#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 10:05:13 2023

@author: e158401a
"""
import cobra
import pandas as pd
import pickle
import diet
from mocbapy.EcosystemModel import create_model, bensolve_default_options
import mocbapy.analysis
from sklearn import metrics
import matplotlib.pyplot as plt
import math
from cobra.flux_analysis import flux_variability_analysis

diet_param = "Western_diet" # Western_diet, Protein_diet or open
host_choice = "sIEC"
plot = True #True to see the Pareto fronts, False for the run to be faster 
o2 = "absent" # is O2 allowed to flow with the food (absent or present)
nbpts_sampling = 10000 #Number of samples for the sampling 
semiconstrained = False #When true, restricts the exchange reaction flux of nutriments without data to -1

def pareto(b, diet_param, host, bacteria_id, enteromet, semiconstrained = semiconstrained, plot=True):
    bacterimet = []
    for met in b.metabolites:
        bacterimet.append((met.id, b))
    all_metabolites = list(set(enteromet + bacterimet))
    metabolic_dict = {x:x[0] for x in all_metabolites}
    if diet_param in ["Western_diet", "High_fiber_diet", "Protein_diet"]:
        host_int = diet.medium_2(model = host.copy(), medium_lumen = diet_param, host_choice = host_choice, block = False, mocba=False, semiconstrained = False)[0]
        growth_host = host_int.optimize().objective_value
        host, diet_dict_1 = diet.medium_2(model = host, medium_lumen = diet_param, host_choice = host_choice, block = False, mocba=True, semiconstrained = False)
        b_int = diet.medium_2(model = b.copy(), medium_lumen = diet_param, host_choice = "bacteria", block = False, mocba=False, semiconstrained = False)[0]
        growth_b = b_int.optimize().objective_value
        b, diet_dict_2 = diet.medium_2(model = b, medium_lumen = diet_param, host_choice = "bacteria", block=False, mocba=True, semiconstrained = False)
        diet_dict = {**diet_dict_1, **diet_dict_2}
        if o2 == "present":        
            diet_dict["o2_e"] = (-10, host.reactions.get_by_id("EX_o2_e").upper_bound)
            print(diet_dict["o2_e"])
        elif o2 == "absent":
            diet_dict["o2_e"] = (0, host.reactions.get_by_id("EX_o2_e").upper_bound)
            #host.reactions.get_by_id("EX_o2_e").bounds = (-10,10)
    elif diet_param == "DM":
        host_int = diet.medium_DM(model = host.copy(), host_choice = host_choice, mocba=False, semiconstrained = semiconstrained)[0]
        growth_host = host_int.optimize().objective_value
        print("growth host alone = "+str(growth_host))
        host, diet_dict_1 = diet.medium_DM(model = host, host_choice = host_choice, mocba=True, semiconstrained = semiconstrained)
        b_int = diet.medium_DM(model = b.copy(), host_choice = "bacteria", mocba=False, semiconstrained = semiconstrained)[0]
        growth_b = b_int.optimize().objective_value
        b, diet_dict_2 = diet.medium_DM(model = b, host_choice = "bacteria", mocba=True, semiconstrained = semiconstrained)
        diet_dict = {**diet_dict_1, **diet_dict_2}
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
        category[0] = "1"
        added_growth_host = maxi_host[0]-growth_host
    if maxi_b[1] > growth_b+(0.01 * growth_b): #growth + 1% of itself to avoid approximations. 
        category[1] = "1"
    if tuple(maxi_host) == tuple(maxi_b) and (category[0] == "1" and category[1] == "1"):
        category[2] = "1"
    category=''.join(category)
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
        #plt.title("AUC = "+str(AUC)+"\n ori_b ="+str(growth_b) +" bwe = "+str(max(y)*growth_b)+"\n ori_e =" +str(growth_host)+ " ewb = "+str(max(x)*growth_host)+"\n category =" + category)
        #plt.title("Pareto front of host - LGG metabolic interaction")
        plt.title("Pareto front of host - bacteria metabolic interaction\n"+str(added_growth_host))
        plt.xlabel("enterocyte's maintenance")
        plt.ylabel(bacteria_id+"'s growth")
        plt.plot(xy['x'], xy['y'], '#ff0000', linestyle="-")
        plt.fill_between(xy['x'], xy['y'], color = "#f08c8c30")
        plt.axhline(y = 1, color = '#1155cc', linestyle = '--', linewidth = 1)
        plt.axvline(x = 1, color = '#1155cc', linestyle = '--', linewidth = 1)
        #plt.savefig("/home/e158401a/Documents/Publi/fig_3/pareto_LGG.pdf", dpi=300)
        plt.show()
    return AUC, xy, ecosys, growth_host, growth_b, diet_dict

def pareto_sampling(ecosys, b, bacteria_id, xy, growth_host, growth_b):
    #conversion from mocbapy to cobra model:
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
    fba = model.optimize() #just to get the reactions names for index of sampling
    sampling_dict = {}
    sampling_dict["reactions"] = fba.fluxes.index
    sumdist = 0        
    for p in xy.index[1:]: #I divide my pareto into segments between 2 points. I denormalize it so that my pareto point can become objective constraints
        x1 = xy.loc[p-1, 'x'] * growth_host
        y1 = xy.loc[p-1, 'y'] * growth_b
        x2 = xy.loc[p, 'x'] * growth_host
        y2 = xy.loc[p, 'y'] * growth_b
        dist = math.hypot(x2-x1, y2-y1)
        sumdist = sumdist + dist
        print(x1, x2, y1, y2, dist)
    d1pt = nbpts_sampling/sumdist
    # All of this was just to get total distance. Could probably be more efficient but will work on it later.
    #Then browse the pareto again to do the sampling.
    ecosystem_biomass = []
    for p in xy.index[1:]:
        print(p)
        x1 = xy.loc[p-1, 'x'] * growth_host
        y1 = xy.loc[p-1, 'y'] * growth_b
        x2 = xy.loc[p, 'x'] * growth_host
        y2 = xy.loc[p, 'y'] * growth_b
        dist = math.hypot(x2-x1, y2-y1)
        nbpnts_xy = round(dist*d1pt)
        print(dist)
        print(nbpnts_xy)
        if nbpnts_xy>0:
            dx = (x2-x1)/nbpnts_xy
            dy = (y2-y1)/nbpnts_xy
            for i in range (0, nbpnts_xy):
                pix = x1+dx*i
                #ecosys.lb[obj_e_i] = pix
                #ecosys.ub[obj_e_i] = pix
                piy = y1+dy*i
                #ecosys.lb[obj_b_i] = piy
                #ecosys.ub[obj_b_i] = piy
                model.reactions.get_by_id('biomass_reactionIEC01b:IEC1907').bounds = (pix, pix)
                model.reactions.get_by_id('Growth:'+b.id).bounds = (piy, piy)
                fba = model.optimize() #One sample
                sampling_dict[str(pix)+"_"+str(piy)] = list(fba.fluxes)
                ecosystem_biomass.append(pix+piy)
                #sampling[str(pix)+"_"+str(piy)] = list(fba.fluxes)            
    sampling = pd.DataFrame(data=sampling_dict)
    sampling.index = sampling["reactions"]
    sampling = sampling.T
    sampling = sampling.drop("reactions")
    return sampling, model

      
def correlation(sampling):
    sampling = sampling.astype(float)
    correlation_reactions = sampling.corr(method ='spearman')
    correlation_reactions = correlation_reactions.fillna(0)
    return correlation_reactions

def short_name(old_name, b_suffixe):
    new_name = old_name.replace(":IEC1907", ":e")
    new_name = new_name.replace(b_suffixe, ":b")
    return new_name

def only_met(reac, b_suffixe):
    met = reac.replace("_e:IEC1907", "")
    met = met.replace("_e:"+b_suffixe[3:], "")
    met = met.replace("EX_", "")
    return met

def oppositeSigns(x, y): #returns boolean. 1 if x y opposite sign, 0 otherwhise
    #return (y >= 0) if (x < 0) else (y < 0) 
    if x < 0:
        return (y > 0) # Because I want a net import 
    elif x > 0: # Because I want a net import 
        return (y < 0)
    else:
        return False

def exchanged_mets(b, bacteria_id, sampling, correlation_reactions,diet_dict):
    ecosystem_biomass = []
    #maxind_b = 0
    maxind_e = 0
    maxind_eco = 0
    maxi_eco = 0
    #maxi_b = 0
    maxi_e = 0
    for i in sampling.index:
        biomass_e = float(i.split("_", 1)[0])
        biomass_b = float(i.split("_", 1)[1])
        ecosystem_biomass.append(biomass_e+biomass_b)
        if biomass_e+biomass_b > maxi_eco:
            maxi_eco = biomass_e+biomass_b
            maxind_eco = i
        if biomass_e > maxi_e:
            maxind_e = i
            maxi_e = biomass_e
        #if biomass_b >= maxi_b:
        #    maxi_b = biomass_b
        #    maxind_b = i            
    e_suffixe = "_e:IEC1907"
    b_suffixe = "_e:"+b.id
    bact_exchanges = []
    for ex in b.exchanges : bact_exchanges.append(ex.id+b_suffixe[2:])
    potential_exchange = {}
    done = []
    for reac in bact_exchanges:
        met = only_met(reac, b_suffixe)
        exmet_e = "EX_"+met+e_suffixe
        exmet_b = "EX_"+met+b_suffixe
        if exmet_e in correlation_reactions.index and exmet_b in correlation_reactions.index:
            #if not math.isnan(correlation_reaction.loc[exmet_e,exmet_b]) and sum(sampling[exmet_e])!=0 and sum(sampling[exmet_b])!=0 and abs(correlation_reaction.loc[exmet_e,exmet_b]) >= 0.5:
            if sum(sampling[exmet_e])!=0 and sum(sampling[exmet_b])!=0 and correlation_reactions.loc[exmet_e,exmet_b] <= -0.5:
                if correlation_reactions.loc[exmet_b,'biomass_reactionIEC01b:IEC1907'] > 0.8:
                    exchange = 0 
                    way = ["0","0"]
                    for s in sampling.index:
                        if oppositeSigns(sampling.loc[s, exmet_e], sampling.loc[s, exmet_b]) :
                            exchange = exchange+1
                            if sampling.loc[s, exmet_e] <0 :
                                way[0] = "1"
                            elif sampling.loc[s, exmet_b]<0:
                                way[1] = "1"
                    if exchange >1 and short_name(exmet_e,b_suffixe)+"|||"+short_name(exmet_b, b_suffixe) not in done:
                        print(exmet_e)
                        potential_exchange[short_name(exmet_e,b_suffixe)+"|||"+short_name(exmet_b, b_suffixe)] = [exchange, way]
                        a = sampling[exmet_e]
                        b = sampling[exmet_b]
                        plt.plot(a, "#e06666")#, label='enterocyte')
                        plt.plot(b, "#93c47d")#, label='bacteria')
                        #plt.plot(ecosystem_biomass, "r", label='Ebiomass')
                        plt.axvline(x = maxind_e, color = "r", linestyle=':')
                        #plt.axvline(x = maxind_e, color = "#e06666", linestyle=':')
                        #plt.axvline(x = maxind_b, color = "#93c47d", linestyle=':')
                        if met+"_e" in list(diet_dict.keys()):
                            plt.axhline(y=diet_dict[met+"_e"][0], color='k', linestyle='--')
                        else:
                            plt.axhline(y=0, color='k', linestyle='--')
                        plt.tick_params(
                            axis='x',          # changes apply to the x-axis
                            which='both',      # both major and minor ticks are affected
                            bottom=False,      # ticks along the bottom edge are off
                            top=False,         # ticks along the top edge are off
                            labelbottom=False) # labels along the bottom edge are off
                        plt.title("evolution of "+met+" exchanges on the pareto front\n"+str(correlation_reactions.loc[exmet_b,'biomass_reactionIEC01b:IEC1907']))
                        plt.legend()
                        plt.savefig("/home/e158401a/Documents/mocba/figures/pareto_exchange_time"+met+"_LGG.png")
                        plt.show()
                done.append(short_name(exmet_e, b_suffixe)+"|||"+short_name(exmet_b, b_suffixe))
    return potential_exchange

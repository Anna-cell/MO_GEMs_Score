#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 16:41:04 2023

@author: e158401a
"""

import cobra
import pandas as pd
import pickle
import analysis
import os
import diet

diet_param = "Western_diet" # Western_diet, Protein_diet or open
host_choice = "sIEC"
plot = True #True to see the Pareto fronts, False for the run to be faster 
o2 = "absent" # is O2 allowed to flow with the food (absent or present)
nbpts_sampling = 100 #Number of samples for the sampling 
semiconstrained = False #When true, restricts the exchange reaction flux of nutriments without data to -1

emblpath = "resources/embl_gems/"
epath = "resources/enterocyte_ori.xml" 

#Import enterocyte model
host = cobra.io.read_sbml_model(epath)
host.objective = host.reactions.get_by_id('biomass_reactionIEC01b')
host.solver = 'cplex'

enteromet = []
for met in host.metabolites:
    enteromet.append((met.id, host))

carve_models = pd.read_csv(emblpath+"/model_list.tsv", sep = "\t", header = 0, index_col=1)

with open("resources/carve_in_vmh.pickle", "rb") as fp:   # EMBL models that are from the gut according to VMH DB
    carve_in_vmh = pickle.load(fp)

#Infer the Pareto front between Lactobacillus rhamnosus GG and the enterocyte     
b = cobra.io.read_sbml_model(emblpath+"models/l/lactobacillus/Lactobacillus_rhamnosus_GG_GG_ATCC_53103.xml")
b.solver = "cplex"
bacteria_id = 'Lactobacillus rhamnosus GG'
AUC, xy, ecosys, growth_host, growth_b, diet_dict = analysis.pareto(b, diet_param, host, bacteria_id, enteromet=enteromet)

#sample solutions on the Pareto front
sampling, model = analysis.pareto_sampling(ecosys, b, bacteria_id, xy, growth_host, growth_b)

#select sampling results from exchange reaction
exchanges = []
for reac in sampling.columns:
    if reac[0:3] == "EX_":
        exchanges.append(reac)
exchanges.append("biomass_reactionIEC01b:IEC1907")
exchanges.append("Growth:Lactobacillus_rhamnosus_GG_GG_ATCC_53103_xml")
sampling_exchanges = sampling[exchanges]
del sampling

# Infer correlation for exchange reaction in the Pareto front
correlation_reactions = analysis.correlation(sampling_exchanges)
# Analyse the metabolites exchanged between both organisms
potential_exchange = analysis.exchanged_mets(b, b.id, sampling_exchanges, correlation_reactions, diet_dict)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 11:03:26 2022

@author: e158401a
"""

import pandas as pd
import pickle


with open("resources/namespace_dicts.pickle", "rb") as fp:   # Unpickling
    namespace_dicts = pickle.load(fp)
namespace_dict_met_inv = {v: k for k, v in namespace_dicts[0].items()}

def no_compartment_id(metabolite):
    a = metabolite
    done = False
    if metabolite[-2:] == "_e":
        metabolite = metabolite.replace('_e','')
        if metabolite != a:
            suffixe='_e'
            done = True
    metabolite = metabolite.replace(('(e)'),'')
    if metabolite != a and not done:
        suffixe='(e)'
        done = True
    metabolite = metabolite.replace(('_c'),'')
    if metabolite != a and not done:
        suffixe='_c'
        done = True
    metabolite = metabolite.replace(('[c]'),'')
    if metabolite != a and not done:
        suffixe='[c]'
        done = True
    metabolite = metabolite.replace(('[m]'),'')
    if metabolite != a and not done:
        suffixe='[m]'
        done = True
    metabolite = metabolite.replace(('[x]'),'')
    if metabolite != a and not done:
        suffixe='[x]'
        done = True
    metabolite = metabolite.replace(('[r]'),'')
    if metabolite != a and not done:
        suffixe='[r]'
        done = True
    metabolite = metabolite.replace(('[e]'),'')
    if metabolite != a and not done:
        suffixe='[e]'
        done = True
    metabolite = metabolite.replace(('[n]'),'')
    if metabolite != a and not done:
        suffixe='[n]'
        done = True
    metabolite = metabolite.replace(('[u]'),'')
    if metabolite != a and not done:
        suffixe='[u]'
        done = True
    metabolite = metabolite.replace(('[luC]'),'')
    if metabolite != a and not done:
        suffixe='[luC]'
        done = True
    metabolite = metabolite.replace(('[bpC]'),'')
    if metabolite != a and not done:
        suffixe='[bpC]'
        done = True
    metabolite = metabolite.replace(('[g]'),'')
    if metabolite != a and not done:
        suffixe='[g]'
        done = True
    metabolite = metabolite.replace(('[l]'),'')
    if metabolite != a and not done:
        suffixe='[l]'
        done = True
    if done == False: 
        print("no suffixe found for", metabolite)
        suffixe = ""
    return metabolite, suffixe

def ex_reacs(model):
    ex_reacs_e = []
    ex_reacs_e_id = []
    for reac in model.reactions:
        #I can't make a model.exchange here because the model doesn't consider (e) reactions as exchange reactions
        find = re.findall(r"EX_",reac.id)
        if find != [] :
            ex_reacs_e.append(reac)
            ex_reacs_e_id.append(reac.id)
    return ex_reacs_e, ex_reacs_e_id
        

AAD_adapted = pd.read_csv("resources/AAD_adapted.tsv", sep="\t", index_col = 0)
with open("resources/mediums_lumen_v2.pickle", "rb") as fp:   # Unpickling
    mediums_lumen = pickle.load(fp)

def medium_2(model, medium_lumen, host_choice, medium_blood = AAD_adapted, block = False, mocba=False, semiconstrained = False):
    constrained_medium_dict = {}
    if host_choice in ["sIEC", "colon"]:
        constrained_medium_dict_bl = {}
        reac_host = [reac.id for reac in model.reactions]
        for i in medium_blood.index:
            reaction = medium_blood.loc[i,host_choice]
            if reaction in reac_host:
                met = list(model.reactions.get_by_id(reaction).metabolites.keys())[0].id
                model.reactions.get_by_id(reaction).bounds = (medium_blood.loc[i,'lb'],medium_blood.loc[i,'ub'])
                constrained_medium_dict_bl[met] = (medium_blood.loc[i,'lb'],medium_blood.loc[i,'ub'])
    for reac in model.exchanges:
        met_e, suffixe = no_compartment_id(list(reac.metabolites.keys())[0].id)
        old_bounds = reac.bounds
        if mocba:
            reac.bounds = (-1000, 1000)
        if met_e in list(mediums_lumen.index):
            constrained_medium_dict[met_e+suffixe] = (-mediums_lumen.loc[met_e, medium_lumen], old_bounds[1])
            if not mocba:
                reac.lower_bound = -mediums_lumen.loc[met_e, medium_lumen]
        elif block == True:
            constrained_medium_dict[met_e+suffixe] = (0, reac._upper_bound)
            if not mocba:
                reac.lower_bound = 0
        elif semiconstrained == True and old_bounds[0] < -1:
            constrained_medium_dict[met_e+suffixe] = (-1, reac._upper_bound)
            if not mocba:
                reac.lower_bound = -1
    return model, constrained_medium_dict

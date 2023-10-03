# MO GEMs Score
This repository contains the code to reproduce results from Lambert et. al. : [fill when published]

This is an analysis of the metabolic interaction between gut microbiota bacteria and a small intestinal epithelial cell using metabolic modeling. Specifically, it uses Multi-objective linear programming to describe trade-offs between organisms's growth (bacteria) or maintenance (human cell). This results in the calculation of an interaction score, revealing the competitive, neutral or mutualist nature of an interaction. 

## Requirements

### CPLEX
To reproduce the results depicted in the publication, this code must be run with the solver CPLEX
### mocbapy
Assembling the individuals models into an ecosystem model is made using [mocbapy](https://gitlab.univ-nantes.fr/mbudinich/mocbapy)
### Python packages
 - pandas
 - cobrapy
 - matplotlib
 - sklearn
 - pickle
 - numpy

## Details
EMBL_enterocyte_scores.py : infer the interaction score between 331 gut bacteria and the enterocyte. Return the interaction categories, added maintenance for the enterocyte in ecosystem, added growth for the bacteria in ecosystem, and that in 3 different diets (open diet, protein diet, western diet)

LGG_enterocyte_interaction.py : samples the pareto front of the interaction between Lactobacillus rhaùnosus GG and the enterocyte in western diet. Infers the model's recations crrelations among these optimal solutions, as well as exchanged metabolites. 

5D_ecosystem_MO.py : Creates an ecosystem of the enterocyte with 4 gut bacteria. Infers the ecosystem score. 

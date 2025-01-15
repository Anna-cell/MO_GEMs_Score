# MO GEMs Score
This repository contains the code to reproduce results from Lambert, A., Budinich, M., Mahé, M., Chaffron, S., & Eveillard, D. (2024). Community metabolic modeling of host-microbiota interactions through multi-objective optimization. Iscience, 27(6). (http://doi.org/10.1016/j.isci.2024.110092)

This is an analysis of the metabolic interaction between gut microbiota bacteria and a small intestinal epithelial cell using metabolic modeling. Specifically, it uses multi-objective linear programming to describe trade-offs between organisms' growth (bacteria) and maintenance (human cell). This results in the calculation of an interaction score, revealing the competitive, neutral, or mutualist nature of an interaction. 

## Requirements

### CPLEX
To reproduce the results depicted in the publication, this code must be run with the solver CPLEX
### mocbapy
Assembling the individual models into an ecosystem model is made using [mocbapy](https://gitlab.univ-nantes.fr/mbudinich/mocbapy)
### Python packages
 - pandas
 - cobrapy
 - matplotlib
 - sklearn
 - pickle
 - NumPy

## Details
EMBL_enterocyte_scores.py: infer the interaction score between 331 gut bacteria and the enterocyte. Return the interaction categories, added maintenance for the enterocyte in ecosystem, added growth for the bacteria in ecosystem, and that in 3 different diets (open diet, protein diet, western diet)

LGG_enterocyte_interaction.py samples the Pareto front of the interaction between Lactobacillus rhamnosus GG and the enterocyte in the Western diet. It infers the model's reaction correlations among these optimal solutions and exchanged metabolites. 

5D_ecosystem_MO.py: Creates an ecosystem of the enterocyte with 4 gut bacteria. Infers the ecosystem score. 

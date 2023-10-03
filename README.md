# MO_GEMs_Score
This repository contains the code to reproduce results from Lambert et. al. : [fill when published]

This is an analysis of the metabolic interaction between gut microbiota bacteria and a small intestinal epithelial cell using metabolic modeling. Specifically, it uses Multi-objective linear programming to describe trade-offs between organisms's growth (bacteria) or maintenance (human cell). This results in the calculation of an interaction score, revealing the competitive, neurral or mutualist nature of an interaction. 

## Requirements

### CPLEX
To reproduce the results depicted in the publication, this code must be run with the solver CPLEX, from IBM
### mocbapy
Assembling the individuals models into an ecosystem model is made using [mocbapy](https://gitlab.univ-nantes.fr/mbudinich/mocbapy)
### Python packages
 - pandas
 - cobrapy
 - matplotlib
 - sklearn
 - pickle
 - numpy

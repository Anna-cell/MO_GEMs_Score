#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:36:17 2023

@author: e158401a
"""
import pickle
import pandas as pd
import numpy as np


# My selection of 4 bacteria, with exact strains instead of refseq

#pareto points
with open("/home/e158401a/Documents/16O-cbapy/sol_mofba_minimal_microbiome_empiric_selection.pickle", "rb") as fp:   # Unpickling
    output = pickle.load(fp)
# sampling
with open("/home/e158401a/Documents/16O-cbapy/sampling10000_ecosys_strains.pickle", "rb") as fp:   # Unpickling
    sampling = pickle.load(fp)
sampling = sampling.sample(n=3000)

points = output[1]
ids = output[2]
#ids = ["enterocyte", "A. muciniphila", "B. xylanisolvens", "F. prausnitzii", "R. bromii"]

growth_ini = output[3]

#######################################################################################
###Calculate a score based on the hypervolume of the pareto linked to the 000 point ###
#######################################################################################

#get rid of slutions with only -1 for one organism and 0 for all the others
points = points[np.all(points >= 0, axis=1)] 
#Add a point for the original growth of each organism, where avery other organism is at 0 (maximal growth when alone)
for i in range(len(growth_ini)):
    ori = [0]*len(growth_ini)
    ori[i] = growth_ini[i]
    points = np.append(points, [ori], axis=0)
#And the original point (every growth to 0) in order to have a comparable form fr each ecosystem. 
ori = [0]*len(growth_ini)
points = np.append(points, [ori], axis=0)

# Normalize each dimension by its organism's initial growth (maximal growth when alone)
for i in range(len(growth_ini)):
    points[:,i] = points[:,i] / growth_ini[i]

#Measure the hypervolume based on the given polyhedron
from scipy.spatial import ConvexHull
volume = ConvexHull(points=points, qhull_options='FA').volume

# Explore the solutions being extreme points of the pareto. Trying to understand which organism influence the score
#the most with a PCA using objective value of the organisms as variables
#Rebuild my data because I must not be normalized by maximal growth when alone here. I also don't need origin points

points = output[1]
ids = ["enterocyte", "A. muciniphila", "B. xylanisolvens", "F. prausnitzii", "R. bromii"]

df_solutions = pd.DataFrame(points, columns=ids)
df_solutions = df_solutions[(df_solutions.iloc[:, :] >= 0).all(axis=1)] #Get rid of -1
df_solutions = df_solutions.reset_index(drop=True)
color = len(df_solutions.index)*["Extreme solution (on the pareto front)"] #Information for PCA : points on pareto

#Solutions where no organism has a null objective
subset = []
for row in df_solutions.index:
    if (df_solutions.loc[row,:] != 0).all():
        subset.append(row)
   
#Identification corresponding to spider plot
color[10] = "Solution 1"
color[13] = "Solution 2"
color[17] = "Solution 3"
color[24] = "Solution 4"


growthes = ['biomass_reactionIEC01b:IEC1907']+["Growth:"+b_id for b_id in output[2][1:]]
sampling = sampling[growthes] #select only biomasses in the sampling
sampling.columns = ids
color = color + len(sampling.index)*["Internal solution"]
df_solutions = pd.concat([df_solutions, sampling])

Ecosystem_biomass = list(df_solutions.transpose().sum())

#pca = PCA(n_components= 3)
#coordinates = pca.fit_transform(df_solutions) #New coordinates (dimension = n_components) for each point ased on pca transformation.

from sklearn.decomposition import PCA
pca = PCA(n_components= 5)
from sklearn.preprocessing import StandardScaler
x_scaled = StandardScaler().fit_transform(df_solutions)
coordinates = pca.fit_transform(x_scaled)

#scree plot
import matplotlib.pyplot as plt
print(pca.explained_variance_ratio_)
PC_values = np.arange(pca.n_components_) + 1
plt.plot(PC_values, pca.explained_variance_ratio_, 'o-', linewidth=2, color='blue')
plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Variance Explained')
plt.savefig("/home/e158401a/Documents/Publi/Fig 4/Fig_4_PCA_5D_screeplot.pdf")
plt.show()

#Find kneepoint to justify only analysing the first components
from kneed import KneeLocator
kn = KneeLocator(PC_values, pca.explained_variance_ratio_, curve='convex', direction='decreasing')
print(kn.knee)
kn.plot_knee_normalized()

pca_df = pd.DataFrame(
    data=coordinates[:,0:2], 
    columns=['PC1', 'PC2'])

import matplotlib.pyplot as plt 
import seaborn as sns

pca_df["Ecosystem biomass"] = Ecosystem_biomass 

sns.scatterplot(
    x='PC1' , 
    y='PC2', 
    hue = "Ecosystem biomass",
    palette= 'YlGn',
    data=pca_df,
    legend=True
    )
plt.title('2D PCA Graph')
#plt.savefig("/home/e158401a/Documents/Publi/Fig 4/Fig_4_PCA_5D_greencolor.pdf")
plt.show()


#loading plot
loadings = pca.components_ #correlation of each component (row) to each feature (col)
n_features = pca.n_features_
#n_components = pca.n_components
feature_names = ids
pc_list = [f'PC{i}' for i in list(range(1, n_features + 1))]
pc_loadings = dict(zip(pc_list, loadings))

loadings_df = pd.DataFrame.from_dict(pc_loadings)
loadings_df['feature_names'] = feature_names
loadings_df = loadings_df.set_index('feature_names')
loadings_df

# 2D
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
sns.set()
 
xs = loadings[0]
ys = loadings[1]

sns.set(rc={'figure.figsize':(8,8)})
sns.set_style(style='white')
sns.scatterplot(
    x='PC1', 
    y='PC2', 
    size = "Ecosystem biomass",#[B*5 for B in Ecosystem_biomass],
    sizes=(20, 100),
    #palette= 'YlGn',
    data=pca_df,
    legend=False,
    hue=color,
    palette = ["black", "grey", "#F37735", "#00B159", "#00AEDB", "#FFC425"],
    hue_order = ["Extreme solution (on the pareto front)", "Internal solution", "Solution 1", "Solution 2", "Solution 3", "Solution 4" ]
    )

a = 3
for i, varnames in enumerate(feature_names):
    plt.scatter(a*xs[i], a*ys[i], s=1)
    plt.arrow(
        0, 0, # coordinates of arrow base
        a*xs[i], # length of the arrow along x
        a*ys[i], # length of the arrow along y
        color='black', 
        head_width=0.01
        )
    plt.text(a*xs[i],a*ys[i], varnames)

plt.xlabel('PC1 '+str(round(pca.explained_variance_ratio_[0], 2)))
plt.ylabel('PC2 '+str(round(pca.explained_variance_ratio_[1], 2)))
plt.title('PCA ecosystem model set of solutions')
#plt.savefig("/home/e158401a/Documents/Publi/Fig 4/Fig_4_PCA_5D.pdf")
plt.show()

import pandas as pd
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('default')
 
pca = PCA(n_components=3)
 
# Fit and transform data
pca_features = pca.fit_transform(x_scaled)
 
# Create dataframe
pca_df = pd.DataFrame(
    data=pca_features, 
    columns=['PC1', 'PC2', 'PC3'])
 
# Initialize the 3D graph
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
 
# Define scaled features as arrays
xdata = pca_df['PC1']
ydata = pca_df['PC2']
zdata = pca_df['PC3']
 
# Plot 3D scatterplot of PCA
ax.scatter3D(
    xdata, 
    ydata, 
    zdata, 
    c=zdata,
    cmap='YlGn', 
    
    alpha=1)
 
# Define the x, y, z variables
loadings = pca.components_
xs = loadings[0]
ys = loadings[1]
zs = loadings[2]
 
# Plot the loadings
for i, varnames in enumerate(feature_names):
    ax.scatter(xs[i], ys[i], zs[i], s=100)
    ax.text(
        xs[i] + 0.1, 
        ys[i] + 0.1, 
        zs[i] + 0.1, 
        varnames)
 
# Plot the arrows
x_arr = np.zeros(len(loadings[0]))
y_arr = z_arr = x_arr
ax.quiver(x_arr, y_arr, z_arr, xs, ys, zs)
 
# Plot title of graph
plt.title('3D Biplot')
 
# Plot x, y, z labels
ax.set_xlabel('Principal component 1', rotation=150)
ax.set_ylabel('Principal component 2')
ax.set_zlabel('Principal component 3', rotation=60)

#ax.view_init(20,0)
plt.show()


growth_alone_eco = pd.DataFrame({"ids" : output[2], "growth_ini" : output[3]})#, "growth_sol" : output[1][1]})
for i in range(len(output[0])):
    if  output[0][i]:
        growth_alone_eco[i] =  output[1][i]

growth_alone_eco.index = growth_alone_eco["ids"]
growth_alone_eco.drop("ids", inplace=True, axis = 1)
growth_alone_eco_T = growth_alone_eco.transpose()
growth_alone_eco_T["Ecosystem_biomass"] = list(growth_alone_eco.sum())

subset = []
for col in growth_alone_eco.columns:
    if (growth_alone_eco[col] != 0).all():
        subset.append(col)

growth_alone_eco_non0 = growth_alone_eco[subset]
growth_alone_eco_T_non0 = growth_alone_eco_non0.transpose()
growth_alone_eco_T_non0["Ecosystem_biomass"] = list(growth_alone_eco_non0.sum())
growth_alone_eco_Ebiomass = growth_alone_eco.transpose()
growth_alone_eco_Ebiomass["Ecosystem_biomass"] = list(growth_alone_eco.sum())


print(growth_alone_eco.sum())
gsum = growth_alone_eco.sum()[1:]

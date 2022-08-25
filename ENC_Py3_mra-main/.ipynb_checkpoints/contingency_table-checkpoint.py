"""
Generates the contingency table for 2 input maps.
"""

import pandas as pd
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
import math

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\Quint\\Documents\\EPA\\Thesis\\ENC_Py3_mra-main\\enrichment_final\\ENC_Py3_mra-main\\"
# Set the case study
case_study = "mmr_2"
# Set the paths to the directories and relevant data
data_path = base_path + "EU_data\\" + case_study + "\\"
output_path = base_path + "EU_output\\" + case_study + "\\"
# Specify the original map (data at time slice 1).
omap_path = data_path + "mmr_landuse_1996_final.tif"
# Specify the actual map (data at time slice 2).
amap_path = data_path + "mmr_landuse_2016_final.tif"
# Specify the masking map.
mask_path = data_path + "mmr_boundary_raster_final.tif"
# Set the land-use class names.
luc_names = ["Agriculture", "Forest", "Industry",
             "Nature", "Recreation", "Urban-formal",
             "Urban-informal", "Urban-villages", "Transportation", "Airport and Harbour",
             "Water","Outside"]

luc = len(luc_names)
# Read in the map for time slice 1.
omap = read_map(omap_path)
# Read in the map for time slice 2.
amap = read_map(amap_path)
# Read in the masking map.
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)

rows = map_dimensions[0]
cols = map_dimensions[1]

def contingency_table(omap, amap, mask, luc, rows, cols):
    cont_table = np.zeros(shape=(luc + 1, luc + 1))
    for i in range(0, rows):
        for j in range(0, cols):
            # Check to make sure the cell is not masked out
            if mask[i, j] != 0:
                # If not masked out, log the presence of the
                # cell between maps 1 & 2.
                x = omap[i, j]
                y = amap[i, j]
                cont_table[x, y] = cont_table[x, y] + 1
    # Determine the total number of each land-use class in map 1.
    for i in range(0, luc):
        for j in range(0, luc):
            x = cont_table[i, j]
            cont_table[i, luc] = cont_table[i, luc] + x
    # Determine the total number of each land-use class in map 2.
    for i in range(0, luc):
        for j in range(0, luc):
            x = cont_table[j, i]
            cont_table[luc, i] = cont_table[luc, i] + x
    # Determine the total number of cells in each map
    for i in range(0, luc):
        x = cont_table[i, luc]
        cont_table[luc, luc] = cont_table[luc, luc] + x
    return cont_table

cont_table = contingency_table(omap, amap, mask, luc, rows, cols)
pd.DataFrame(cont_table).to_csv("cont_table.csv")
#cont_table.to_csv(base_path+'cont_table.csv')

np.savetxt("cont_table.csv", cont_table, delimiter=",")

file = open(base_path+"cont_table.txt", "w+")
#Saving the array in a text file
content = str(cont_table)
file.write(content)
file.close()
 
print(cont_table)
np.set_printoptions(threshold = False)
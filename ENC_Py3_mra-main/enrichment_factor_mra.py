from enrichment_factor import ef
# Modules
# Allows for multi-dimensional array handling.
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)

# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
import math

base_path = "C:\\Users\\Gamelab\\Desktop\\RT\\Others\\Thesis\\Thesis_coding\\ABM\\ENC_Py3_mra\\"
# Set the case study
case_study = "Amsterdam"
# Set the paths to the directories and relevant data
data_path = base_path + "EU_data\\" + case_study + "\\"
output_path = base_path + "EU_output\\" + case_study + "\\"
# Specify the original map (data at time slice 1).
omap_path = r'enrichment_final/ENC_Py3_mra-main/mmr_landuse_1996_original.tif'
# Specify the actual map (data at time slice 2).
amap_path = r'enrichment_final/ENC_Py3_mra-main/mmr_landuse_2016_original.tif'
# Specify the masking map.
mask_path = r'enrichment_final/ENC_Py3_mra-main/mmr_boundary_raster_final.tif'
# Set the land-use class names.
luc_names = ["Outside","Agriculture", "Forest", "Industry",
             "Nature", "Recreation", "Urban-formal",
             "Urban-informal", "Urban-villages", "Transportation", "Airport and Harbour",
             "Water"]

# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 3
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 8

# Set the minimum required percentage of conversions.
min_convo_rate = 0.025
# Set the minimum required Enrichment Factor value at distance 1.
min_EF_1 = 0.00
# Set the significance limit (recommended value is 1.96
z_limit = 1.96

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

# Determine the distances that will be analysed using the module considered
# distances.
temp = considered_distances(max_distance)
# Store the list of considered distances as a variable.
cd = temp[0]
# Store the total number of distances considered
cdl = temp[1]
# Determine the maximum neighbourhood size (unit) from considered distances
N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
N = []
for c in range(0, max_distance):
    N.append(N_all[c])




enriched_factors = ef(luc, max_d=8, cdl=cdl, cd=cd, N=N, omap=omap, amap=amap, mask=mask, row=rows, col=cols)

enriched_factors.to_csv('enrichment_MMR_5.csv')

np.set_printoptions(threshold = False)


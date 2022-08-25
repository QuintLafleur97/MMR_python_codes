"""
Stage 2 of the Empirical Neighbourhood Calibration method:
Categorise the different interactions per group based on the
specified data.
"""

# Modules
# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Generate a contingency table for the data.
from contingency_table import contingency_table
# Calculate the enrichment factor.
from enrichment_factor_Q_2 import ef
# Log scale the enrichment factor values (base 10).
from log_scale_ef import log_scale_ef

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
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 4
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 8


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

# Count the presence of each land-use class in the actual map. This is
# used in the calculation of area-weighted average clumpiness across the
# active classes.
luc_count = [0] * luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count[amap[i, j]] = luc_count[amap[i, j]] + 1

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
points_path=output_path+"Rules\\"+"points.txt"
tails_path=output_path+"Rules\\"+"tails.txt"
# Read in the points file of inertial and conversion points.
points = np.loadtxt(points_path)
# Read in the tails file of inertial and conversion tails.
tails = np.loadtxt(tails_path)

# Data analysis --------------------------------------------------------------
# Generate a contingency table for the data.
cont_table = contingency_table(omap, amap, mask, luc, rows, cols)
# Determine the rates of inertia and conversion.
ic_rates = np.zeros(shape=(luc, luc))
for i in range(0, luc):
    for j in range(0, luc):
        if i == j:
            if cont_table[i, luc] > 0:
                ic_rates[i, j] = cont_table[i, j] / cont_table[i, luc]
        else:
            conversions = abs(float(cont_table[j, j]) - float(cont_table[luc, j]))
            if conversions > 0:
                ic_rates[i, j] = float(cont_table[i, j]) / float(conversions)
# Determine the Enrichment Factor values.

data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)
print(data_ef)
log_data_ef = log_scale_ef(data_ef, 10, luc, act, pas, max_distance)

# Categorisation --------------------------------------------------------------
# This is implemented for 3 categories per interaction type (low = 1, med = 2,
# high = 3).

# Specify the category levels.
# The inertia rate levels for setting inertial points.
high_inertia_rate = 0.95
mid_inertia_rate = 0.90
# The conversion rate levels for setting conversion points.
high_conversion_rate = 0.5
mid_conversion_rate = 0.1
# The enrichment factor levels for setting the tail values.
high_ef = 1.0
mid_ef = 0.5
# Generate two metrics to track the levels for the points and tails.
point_levels = np.zeros(shape=(luc,act))
tail_levels = np.zeros(shape=(luc,act))

# Analyse the data, set the levels.
for i in range(0, act):
    for j in range(0, luc):
        # Analyse the points included for calibration
        if points[j, i] == 1:
            # Evaluate inertia points.
            if i + pas == j:
                point_levels[j, i] = 1
                inertia_rate = ic_rates[j, i + pas]
                if cont_table[i + pas, luc] > cont_table[luc, i + pas]:
                    point_levels[j, i] = 1
                elif inertia_rate > high_inertia_rate:
                    point_levels[j, i] = 3
                elif inertia_rate > mid_inertia_rate:
                    point_levels[j, i] = 2
            # Evaluate the conversion points.
            else:
                conversion_rate = ic_rates[j, i + pas]
                if conversion_rate > high_conversion_rate:
                    point_levels[j, i] = 3
                elif conversion_rate > mid_conversion_rate:
                    point_levels[j, i] = 2
                else:
                    point_levels[j, i] = 1
        # Analyse the tails included for calibration.
        if tails[j, i] == 1:
            # Evaluate the inertia tails.
            if i + pas == j:
                tail_levels[j, i] = 1
                if log_data_ef[1, j, i] > high_ef:
                    tail_levels[j, i] = 3
                elif log_data_ef[1, j, i] > mid_ef:
                    tail_levels[j, i] = 2
            # Evaluate the conversion tails.
            else:
                tail_levels[j, i] = 1
                if log_data_ef[1, j, i] > high_ef:
                    tail_levels[j, i] = 3
                elif log_data_ef[1, j, i] > mid_ef:
                    tail_levels[j, i] = 2

# Save the output.
tail_levels_file = output_path + "Rules\\tail_levels2.txt"
point_levels_file = output_path + "Rules\\point_levels2.txt"
np.savetxt(tail_levels_file, tail_levels, fmt="%d")
np.savetxt(point_levels_file, point_levels, fmt="%d")

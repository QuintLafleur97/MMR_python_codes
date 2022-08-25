"""
This module can be used to calculate both Cohen's Kappa and Kappa Simulation,
a derivative that focuses on the crisp agreement of transitions relative to a
baseline of random allocation (Van Vliet et al. 2011).
"""


import numpy as np


def kappa(map1, map2, mask):
    # Determine the map dimensions and number of land-use classes.
    shape_map1 = np.shape(map1)
    row = shape_map1[0]
    column = shape_map1[1]
    luc = np.amax(map1) + 1
    # Determine the total number of cells to be considered.
    total_cells = 0
    for i in range(0, row):
        for j in range(0, column):
            x = mask[i, j]
            if x != 0:
                total_cells = total_cells + 1
    # Initialise an array to store the observed agreement probability.
    Po_array = np.zeros(shape=luc)
    # Initialise a set of arrays to store the expected agreement probability,
    # for both maps, and then combined.
    Pe1_array = np.zeros(shape=luc)
    Pe2_array = np.zeros(shape=luc)
    Pe_array = np.zeros(shape=luc)
    # Initialise an array to store the maximum possible agreement probability.
    Pmax_array=np.zeros(shape=luc)
    # Analyse the agreement between the two maps.
    for i in range(0, row):
        for j in range(0, column):
            if mask[i,j] != 0:
                x = map1[i, j]
                y = map2[i, j]
                # Amend the expected array count
                Pe1_array[x]=Pe1_array[x]+1
                Pe2_array[y]=Pe2_array[y]+1
                # If there is agreement, amend the observed count
                if x == y:
                        Po_array[x] = Po_array[x] + 1
    # Convert to probabilities.
    Po_array[:] = [x/total_cells for x in Po_array]
    Pe1_array[:] = [x/total_cells for x in Pe1_array]
    Pe2_array[:] = [x/total_cells for x in Pe2_array]
    # Now process the arrays to determine the maximum and expected
    # probabilities.
    for i in range(0, luc):
        Pmax_array[i] = min(Pe1_array[i], Pe2_array[i])
        Pe_array[i] = Pe1_array[i]*Pe2_array[i]
    # Calculate the values of probability observed, expected, and max.
    Po = np.sum(Po_array)
    Pmax = np.sum(Pmax_array)
    Pe = np.sum(Pe_array)
    # Now calculate the Kappa histogram and Kappa location.
    Khist = (Pmax - Pe)/(1 - Pe)
    Kloc = (Po - Pe)/(Pmax - Pe)
    # Finally, calculate Kappa.
    Kappa=Khist*Kloc
    # Return the value of Kappa.
    return Kappa


def ksim(omap, map1, map2, mask):
    # Determine the map dimensions and number of land-use classes.
    shape = np.shape(map1)
    row = shape[0]
    column = shape[1]
    luc = np.amax(map1) + 1
    # Determine the total number of cells to be considered.
    total_cells = 0
    for i in range(0, row):
        for j in range(0, column):
            x = mask[i, j]
            if x != 0:
                total_cells = total_cells + 1
    # Initialise an array to store the observed and expected agreement for the
    # transitions between the two maps.
    Po_o_array = np.zeros(shape=luc)
    Po_a_array = np.zeros(shape=luc)
    Pe1trans_array = np.zeros(shape=luc**2)
    Pe2trans_array = np.zeros(shape=luc**2)
    Petrans_array = np.zeros(shape=luc)
    Pmaxtrans_array = np.zeros(shape=luc)
    # Evaluate the maps via couting.
    for i in range(0, row):
        for j in range(0, column):
            if mask[i, j] != 0:
                x=omap[i, j]
                y=map1[i, j]
                z=map2[i, j]
                
                Po_o_array[x] = Po_o_array[x] + 1
                
                Pe1trans_array[x*luc + y] = Pe1trans_array[x*luc + y] + 1
                Pe2trans_array[x*luc + z] = Pe2trans_array[x*luc + z] + 1
                
                if y==z:
                    Po_a_array[z] = Po_a_array[z] + 1
    # Convert the counts to proportions
    Po_o_array[:] = [x/total_cells for x in Po_o_array]
    Po_a_array[:] = [x/total_cells for x in Po_a_array]
    Pe1trans_array[:] = [x/total_cells for x in Pe1trans_array]
    Pe2trans_array[:] = [x/total_cells for x in Pe2trans_array]
    # Analyse the observed agreements for calculation of the expected
    # agreements.
    for i in range(0, luc):
        for j in range(0, luc):
            if Po_o_array[i] > 0:
                Pe1trans_array[i*luc + j] = (
                    Pe1trans_array[i * luc + j] / Po_o_array[i]
                )
                Pe2trans_array[i*luc + j] = (
                    Pe2trans_array[i * luc + j] / Po_o_array[i]
                )

    for i in range(0, luc):
        for j in range(0, luc):
            Petrans_array[i] = (
                Petrans_array[i] + Po_o_array[i] *
                Pe1trans_array[luc * i + j] * Pe2trans_array[luc * i + j]
            )
            Pmaxtrans_array[i] = (
                Pmaxtrans_array[i] + Po_o_array[i] *
                min(Pe1trans_array[luc * i + j], Pe2trans_array[luc * i + j])
            )
    # Calculate the values of observed expected, and max transition agreement.
    Po_a = np.sum(Po_a_array)
    Petrans = np.sum(Petrans_array)
    Pmaxtrans = np.sum(Pmaxtrans_array)
    # Print error feedback if the solution cannot be evaluated.
    if Petrans == 1:
        return 'Calculation results in undefined solution'
        
    if Petrans == Pmaxtrans:
        return 'Calculation results in undefined solution'
    # Calculate the histogram and location of Kappa Simulation.
    Ktransition = (Pmaxtrans - Petrans)/(1 - Petrans)
    Ktransloc = (Po_a - Petrans)/(Pmaxtrans - Petrans)
    # Finally, calculate Kappa Simulation, and return the value.
    KSIM = Ktransition*Ktransloc
    return KSIM


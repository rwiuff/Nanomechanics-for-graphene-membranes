# =====================================================================
#
# By Tue Gunst & Rasmus Wiuff
#
# Script for correcting constrain tags to ommit ATKPython bug
#
# =====================================================================

# -------------------------------------------------------------
# Import Libraries
# -------------------------------------------------------------
from NanoLanguage import *
import numpy as np

# -------------------------------------------------------------
# Set layermode
# -------------------------------------------------------------
lm = 1
#lm = 0

# -------------------------------------------------------------
# Load File
# -------------------------------------------------------------
path = 'Sheet.hdf5'
bulk_configuration = nlread(path, BulkConfiguration)[-1]

# -------------------------------------------------------------
# Create Substrate tag
# -------------------------------------------------------------
if lm == 1:
    # Get atom coordinates
    coor = bulk_configuration.fractionalCoordinates()
    sub = np.array([])
    # Sort coordinates for 1 layer
    for i in range(coor.shape[0]):
        if coor[i, 2] < 0.50:
            sub = np.append([sub], [i])
    # Convert indices to list and define tag with constraints
    sub = sub.astype(int)
    sub = sub.tolist()
    bulk_configuration.addTags('Substrate', sub)

    constraints = bulk_configuration.indicesFromTags('Substrate')
elif lm == 0:
    # Define tag with constraints
    constraints = bulk_configuration.indicesFromTags('Substrate')


# -------------------------------------------------------------
# Fix Constraints
# -------------------------------------------------------------
# Make indices of atoms in configuration and get constraints.
atomindices = numpy.arange(bulk_configuration.numberOfAtoms())
# Move constrained atoms to end of list.
atomindices_new = numpy.delete(atomindices, constraints)
atomindices_new = numpy.append(atomindices_new, constraints)
constraints_new = atomindices[len(atomindices) - len(constraints)::]
new_xyz = bulk_configuration.cartesianCoordinates().inUnitsOf(Angstrom)
new_xyz = new_xyz[atomindices_new]
# Setup new configuration.
bulk_configuration._changeAtoms(positions=new_xyz * Angstrom)
constraints = constraints_new
bulk_configuration.addTags('FixedFinal', constraints)

# -------------------------------------------------------------
# Make Layer tags
# -------------------------------------------------------------
if lm == 1:
    coor = bulk_configuration.fractionalCoordinates()
    l1 = np.array([])
    l2 = np.array([])
    for i in range(coor.shape[0]):
        if coor[i, 2] < 0.50:
            l2 = np.append([l2], [i])
        elif coor[i, 2] > 0.50:
            l1 = np.append([l1], [i])
    l1 = l1.astype(int)
    l2 = l2.astype(int)
    l1 = l1.tolist()
    l2 = l2.tolist()
    bulk_configuration.addTags('Layer1', l1)
    bulk_configuration.addTags('Layer2', l2)
# Save configuration with fixed tags
nlsave('SheetFixed.hdf5', bulk_configuration)
nlprint(bulk_configuration)

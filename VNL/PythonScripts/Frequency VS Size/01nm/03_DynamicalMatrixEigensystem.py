# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Import Libraries
# -------------------------------------------------------------
#import pandas as pd
import numpy

# -------------------------------------------------------------
# Import DynamicalMatrix
# -------------------------------------------------------------
path = 'DynamicalMatrix.hdf5'
configuration = nlread(path, BulkConfiguration)[-1]
dynamical_matrix = nlread(path, DynamicalMatrix)[-1]

# -------------------------------------------------------------
# Extract phononEigensystem
# -------------------------------------------------------------
# Get dimensions
n_modes = (len(configuration) - len(dynamical_matrix.constraints())) * 3

# Find constraints.
if len(dynamical_matrix.constraints()) > 0:
    constrained_displacements = numpy.sort(list(3 * dynamical_matrix.constraints()) + list(
        3 * dynamical_matrix.constraints() + 1) + list(3 * dynamical_matrix.constraints() + 2)).astype(int)

# Calculate eigevalues and eigenvectors from dynamical matrix.
qpoint = [0, 0, 0]
frequency_list = numpy.zeros((n_modes), dtype=float)
eigenvalues, eigenvectors = dynamical_matrix.phononEigensystem(qpoint)
frequency_list[:] = eigenvalues.inUnitsOf(eV)
for mode in numpy.arange(n_modes):
    # Remove constrained displacements.
    if len(dynamical_matrix.constraints()) > 0:
        eigenvector = numpy.delete(
            eigenvectors[:, mode], constrained_displacements)
    else:
        eigenvector = eigenvectors[:, mode]
numpy.savez_compressed(
    'EigenSystem', EigenValues=eigenvalues, EigenVectors=eigenvectors)
loaded = numpy.load('EigenSystem.npz')
print(numpy.array_equal(eigenvalues, loaded['EigenValues']))
print(numpy.array_equal(eigenvectors, loaded['EigenVectors']))

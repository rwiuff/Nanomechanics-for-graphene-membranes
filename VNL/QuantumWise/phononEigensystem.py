# --------------------------------------------------------------------------- #
# Load libraries
from NanoLanguage import *
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Choose file
myfile = 'DynamicalMatrix.hdf5'
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Load configuration
configuration = nlread(myfile, BulkConfiguration)[-1]
# Load DynamicalMatrix
dynamical_matrix = nlread(myfile, DynamicalMatrix)[-1]
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Calculate eigevalues and eigenvectors from dynamical matrix at first q-point
eigenvalues, eigenvectors = dynamical_matrix.phononEigensystem([0, 0, 0])
# --------------------------------------------------------------------------- #

print("Number of eigenvalues: {}".format(eigenvalues.shape))
print("Shape of eigenvector matrix: {}".format(eigenvectors.shape))

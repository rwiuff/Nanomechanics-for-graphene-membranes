# --------------------------------------------------------------------------- #
# Load libraries
from NanoLanguage import *
from pylab import *
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Choose file
myfile = 'DynamicalMatrix.hdf5'
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Load configuration with calculator
configuration = nlread(myfile, BulkConfiguration)[-1]
# Load DynamicalMatrix
dynamical_matrix = nlread(myfile, DynamicalMatrix)[-1]
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Calculate eigevalues and eigenvectors from dynamical matrix.
eigenvalues, eigenvectors = dynamical_matrix.phononEigensystem([0,0,0])
# --------------------------------------------------------------------------- #

print(eigenvalues.shape)
print(eigenvectors.shape)

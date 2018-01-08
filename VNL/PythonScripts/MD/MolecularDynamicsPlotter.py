# -------------------------------------------------------------
# Load libraries, ProjectedPhononBandsDisplacement and switch
# matplotlib backend
# -------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.lines as mlines
from NanoLanguage import *
from pylab import *
import pickle
import numpy as np
import numpy
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from ConfigDisplacementFunctions import ChangePositionsBulk, myPhononModeDisplacedConfig

# -------------------------------------------------------------
# Load data from dynamical matrices
# -------------------------------------------------------------
ZNumber = 2
if ZNumber >= 10:
    # Choose file
    myfile = '{}DynamicalMatrix.hdf5'.format(ZNumber)
else:
    myfile = '0{}DynamicalMatrix.hdf5'.format(ZNumber)
# Load configuration with calculator
configuration = nlread(myfile, BulkConfiguration)[-1]

# Load DynamicalMatrix
dynamical_matrix = nlread(myfile, DynamicalMatrix)[-1]

# Vibrational state to project onto
n_modes = (len(configuration) - len(dynamical_matrix.constraints())) * 3

# -------------------------------------------------------------
# Load data from clamped reference dynamical matrix
# -------------------------------------------------------------
RFDM = '05nmDynamicalMatrix.hdf5'
RFconfiguration = nlread(RFDM, BulkConfiguration)[-1]
RFdynamical_matrix = nlread(RFDM, DynamicalMatrix)[-1]
RFn_modes = ((len(RFconfiguration) - len(
    RFdynamical_matrix.constraints())) * 3)

# Display loaded matrices and modes
print('+------------------------------------------------------+')
pstring = "| The File {} contains {:4d} modes |".format(myfile, int(n_modes))
print(pstring)
print('+------------------------------------------------------+')

print('+-------------------------------------------------------+')
pstring = "| The File {} contains {:4d} modes |".format(RFDM, int(RFn_modes))
print(pstring)
print('+-------------------------------------------------------+')

# -------------------------------------------------------------
# Calculate the phononEigensystems
# -------------------------------------------------------------
constraints = dynamical_matrix.constraints()
constraints = map(int, constraints)
DMvibrational_frequencies, DMvibrational_states = RFdynamical_matrix.phononEigensystem(constrained_atoms=constraints)

# -------------------------------------------------------------
# Displace the configurations
# -------------------------------------------------------------
inverse_characteristic_length, sigma, configuration_new, dRt = myPhononModeDisplacedConfig(configuration=configuration, vibrational_frequencies = DMvibrational_frequencies, vibrational_states=DMvibrational_states, phonon_mode=0)
nlsave('test.hdf5',configuration_new)
# -------------------------------------------------------------
# Plot the energies
# -------------------------------------------------------------

# -------------------------------------------------------------
# Load libraries, ProjectedPhononBandsDisplacement and switch
# matplotlib backend
# -------------------------------------------------------------
import matplotlib.pyplot as plt
from NanoLanguage import *
from pylab import *
import pickle
import numpy as np

# -------------------------------------------------------------
# Create data arrays
# -------------------------------------------------------------
configuration = np.array([])
vibrationalmodes = np.array([])

# -------------------------------------------------------------
# Load data from clamped reference dynamical matrix
# -------------------------------------------------------------
RFDM = '05nmDynamicalMatrix.hdf5'
RFconfiguration = nlread(RFDM, BulkConfiguration)[-1]
RFdynamical_matrix = nlread(RFDM, DynamicalMatrix)[-1]
nof = 3
myfile = np.array([])
for i in range(nof + 1):
    if i == nof:
        break
    elif i >= 9:
        # Choose file
        myfile = np.append(myfile,
                           '{}SheetVib.hdf5'.format(i + 1))
    else:
        myfile = np.append(myfile,
                           '0{}SheetVib.hdf5'.format(i + 1))

RFSV = '05nmSheetVib.hdf5'
RFvibrationalmodes = nlread(RFSV, VibrationalMode)[-1]

print("+========================+")
print("|    Saving datafiles    |")
print("--------------------------")
with open('ClampedModes.pickle', 'wb') as handle:
    pickle.dump(anti_projection, handle,
                protocol=pickle.HIGHEST_PROTOCOL)
print("| (4/5): ClampedModes |")
with open('01Modes.pickle', 'wb') as handle:
    pickle.dump(qpoints, handle,
                protocol=pickle.HIGHEST_PROTOCOL)
print("| (1/5): 01Modes        |")
with open('02Modes.pickle', 'wb') as handle:
    pickle.dump(frequency_list, handle,
                protocol=pickle.HIGHEST_PROTOCOL)
print("| (2/5): 02Modes  |")
with open('03Modes.pickle', 'wb') as handle:
    pickle.dump(projection, handle,
                protocol=pickle.HIGHEST_PROTOCOL)
print("| (3/5): 03Modes      |")
print("+========================+")
quit()

# -------------------------------------------------------------
# Load libraries
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

ClampedModes = {}
ZModes = {}
T = 300000
Index = [0, 1, 3, 5, 6]
coor = np.array([])
for i in range(5):
    print('+----------------------------------+')
    pstring = "| Extracting clamped mode {:2d} of {:2d} |".format(
        int(1 + i), int('5'))
    print(pstring)
    print('+----------------------------------+')
    VM = RFvibrationalmodes.movie(
        mode_index=Index[i], temperature=T * Kelvin)
    ClampedModes[i] = VM.lastImage()
    VM = np.array([])
RFvibrationalmodes = np.array([])
ModeIndex = [2, 3, 5, 7, 8]
ImageIndex = np.array([[9, 19, 9, 9, 19],
                       [9, 9, 19, 19, 9],
                       [19, 19, 19, 19, 9]])
coor = np.array([])
c = 0
for j in range(3):
    vibrationalmodes = nlread(myfile[j], VibrationalMode)[-1]
    for i in range(5):
        c = c + 1
        print('+----------------------------------+')
        pstring = "| Extracting mode         {:2d} of {:2d} |".format(
            int(c), int('15'))
        print(pstring)
        print('+----------------------------------+')
		# Extract the trajectory from VibrationalMode at mode i and
		# temperature T
        VM = vibrationalmodes.movie(
            mode_index=ModeIndex[i], temperature=T * Kelvin)
		# Extract configuration from image j,i in trajectory
        ZModes[c] = VM.image(image_index=ImageIndex[j, i])
        VM = np.array([])


print("+======================+")
print("|   Saving datafiles   |")
print("------------------------")
with open('ClampedModes.pickle', 'wb') as handle:
    pickle.dump(ClampedModes, handle,
                protocol=pickle.HIGHEST_PROTOCOL)
print("| (1/2): ClampedModes  |")
with open('ZModes.pickle', 'wb') as handle:
    pickle.dump(ZModes, handle,
                protocol=pickle.HIGHEST_PROTOCOL)
print("| (2/2): ZModes        |")
print("+======================+")

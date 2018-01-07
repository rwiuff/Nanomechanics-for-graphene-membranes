# -------------------------------------------------------------
# Load libraries, ProjectedPhononBandsDisplacement and switch
# matplotlib backend
# -------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib import gridspec
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

# -------------------------------------------------------------
# Define colormaps
# -------------------------------------------------------------
cmap, norm = cm.get_cmap('brg'), None  # 'hot','brg','seismic'

# -------------------------------------------------------------
# Plot RMS
# -------------------------------------------------------------
fig = plt.figure(figsize=(5, 4))
gs = gridspec.GridSpec(4, 5)
gs.update(wspace=0.03, hspace=0.03)
ax1 = {}
ax2 = {}
Index = [0, 1, 3, 5, 6]
coor = np.array([])
for i in range(5):
    print('+--------------------------------+')
    pstring = "| Plotting clamped plot {:2d} of {:2d} |".format(
        int(1 + i), int('5'))
    print(pstring)
    print('+--------------------------------+')
    ax1[i] = plt.subplot(gs[0, i])
    VM = RFvibrationalmodes.movie(
        mode_index=Index[i], temperature=300000 * Kelvin)
    coor = VM.image(19)
    coor = coor.cartesianCoordinates().inUnitsOf(Ang)
    VM = np.array([])
    x = coor[:, 0]
    y = coor[:, 1]
    z = coor[:, 2]
    z = (coor[:, 2]) - 10
    plt.tricontourf(x, y, z, 200)
    plt.ylim(-31.1505, 31.1505)
    plt.xlim(50.464, 112.765)

for i in range(1, 5):
    plt.setp(ax1[i].get_yticklabels(), visible=False)
    plt.subplots_adjust(left=0.07, bottom=None, right=0.88, top=None,
                        wspace=None, hspace=None)
for i in range(5):
    plt.setp(ax1[i].get_xticklabels(), visible=False)
    plt.subplots_adjust(left=0.07, bottom=None, right=0.88, top=None,
                        wspace=None, hspace=None)

RFvibrationalmodes = np.array([])
Index = [2, 3, 5, 7, 8]
coor = np.array([])
c = 0
for j in range(3):
    vibrationalmodes = nlread(myfile[j], VibrationalMode)[-1]
    for i in range(5):
        c = c + 1
        print('+--------------------------------+')
        pstring = "| Plotting mode plot    {:2d} of {:2d} |".format(
            int(c), int('15'))
        print(pstring)
        print('+--------------------------------+')
        ax2[c - 1] = plt.subplot(gs[j + 1, i])
        VM = vibrationalmodes.movie(
            mode_index=Index[i], temperature=300000 * Kelvin)
        coor = VM.image(19)
        coor = coor.cartesianCoordinates().inUnitsOf(Ang)
        rmindex = np.array([])
        for i in range(coor.shape[0]):
            if coor[i, 2] < 10:
                rmindex = np.append(rmindex, [i])
        coor = np.delete(coor, rmindex, axis=0)
        VM = np.array([])
        x = coor[:, 0]
        y = coor[:, 1]
        z = (coor[:, 2]) - 11.8747
        plt.tricontourf(x, y, z, 200)
        plt.ylim(-30.5601, 30.5601)
        plt.xlim(20.5441, 82.8171)

rangeindex = [0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 10, 0, 0, 0, 0]
for i in range(15):
    if i != rangeindex[i]:
        plt.setp(ax2[i].get_yticklabels(), visible=False)
        plt.subplots_adjust(left=0.07, bottom=None, right=0.88, top=None,
                            wspace=None, hspace=None)
for i in range(10):
    plt.setp(ax2[i].get_xticklabels(), visible=False)
    plt.subplots_adjust(left=0.07, bottom=None, right=0.88, top=None,
                        wspace=None, hspace=None)


cax = plt.axes([0.89, 0.1, 0.03, 0.8])
cb = plt.colorbar(cax=cax)
cb.set_label('Displacement [$\AA$]', fontsize=12)
plt.show()
quit()
# -------------------------------------------------------------
# Show or save plots
# -------------------------------------------------------------
showsave = 0
menu = np.array(["Show plot                  ",
                 "Save plot                  "])
while showsave == 0:
    print("+================================+")
    print("|       Show or save plot?       |")
    print("|--------------------------------|")
    for i in range(len(menu)):
        print("| {:d}. {:s} |".format(i + 1, menu[i]))
    print("+================================+")
    while not (np.any(showsave == np.arange(len(menu)) + 1)):
        while True:
            try:
                showsave = float(input("Choose an option: "))
                break
            except ValueError:
                print("Only integers accepted")
                pass
            except NameError:
                print("Only integers accepted")
            except TypeError:
                print("Only integers accepted")
if showsave == 1:
    plt.show()
elif showsave == 2:
    ylim(ymin, ymax)
    print("+=====================================+")
    print("|              Saving plots           |")
    print("|-------------------------------------|")
    print("| (1/2): ZetaModeProjections.eps |")
    savefig('ZetaModeProjections.eps')
    print("| (2/2): ZetaModeProjections.png |")
    savefig('ZetaModeProjections.png')
    print("+=====================================+")

# -------------------------------------------------------------
# Load libraries, ProjectedPhononBandsDisplacement and switch
# matplotlib backend
# -------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as colors
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
# Load data from pickle files
# -------------------------------------------------------------
ZModes = {}
ClampedModes = {}

print("+=======================+")
print("|   Loading datafiles   |")
print("|-----------------------|")
with open('ClampedModes.pickle', 'rb') as handle:
    ClampedModes = pickle.load(handle)
print("| (1/2): Clamped Modes  |")
with open('ZModes.pickle', 'rb') as handle:
    ZModes = pickle.load(handle)
print("| (2/2): ZModes         |")
print("+=======================+")

# -------------------------------------------------------------
# Define colormaps
# -------------------------------------------------------------
cmap = cm.get_cmap('brg')
v_min = 0
v_max = 0
for i in range(len(ClampedModes)):
    if v_min > np.min((
            ClampedModes[i].cartesianCoordinates().inUnitsOf(Ang)[:, 2]) - 10):
        v_min = np.min(
            (ClampedModes[i].cartesianCoordinates().inUnitsOf(Ang)[:, 2]) - 10)
    else:
        v_min = v_min
for i in range(1, len(ZModes)):
    coor = ZModes[i]
    coor = coor.cartesianCoordinates().inUnitsOf(Ang)
    rmindex = np.array([])
    for i in range(coor.shape[0]):
        if coor[i, 2] < 10:
            rmindex = np.append(rmindex, [i])
    coor = np.delete(coor, rmindex, axis=0)
    z = (coor[:, 2]) - 11.8747
    if v_min > np.min(z):
        v_min = np.min(z)
    else:
        v_min = v_min
for i in range(len(ClampedModes)):
    if v_max < np.max((
            ClampedModes[i].cartesianCoordinates().inUnitsOf(Ang)[:, 2]) - 10):
        v_max = np.max(
            (ClampedModes[i].cartesianCoordinates().inUnitsOf(Ang)[:, 2]) - 10)
    else:
        v_max = v_max
for i in range(1, len(ZModes)):
    if v_max < np.max((
            ZModes[i].cartesianCoordinates().inUnitsOf(Ang)[:, 2]) - 11.8747):
        v_max = np.max(
            (ZModes[i].cartesianCoordinates().inUnitsOf(Ang)[:, 2]) - 11.8747)
    else:
        v_max = v_max
levels = np.linspace(v_min, v_max, 800)
# -------------------------------------------------------------
# Plot RMS
# -------------------------------------------------------------
Index = [0, 1, 3, 5, 6]
np.set_printoptions(threshold='nan')
fig = plt.figure(figsize=(10.15, 8.8))
gs = gridspec.GridSpec(4, 5)
gs.update(wspace=0.03, hspace=0.2)
ax1 = {}
ax2 = {}
coor = np.array([])
for i in range(5):
    print('+--------------------------------+')
    pstring = "| Plotting clamped plot {:2d} of {:2d} |".format(
        int(1 + i), int('5'))
    print(pstring)
    print('+--------------------------------+')
    ax1[i] = plt.subplot(gs[0, i])
    coor = ClampedModes[i]
    coor = coor.cartesianCoordinates().inUnitsOf(Ang)
    x = coor[:, 0]
    y = coor[:, 1]
    z = coor[:, 2]
    z = (coor[:, 2]) - 10
    plt.tricontourf(x, y, z, vmin=v_min, vmax=v_max, levels=levels)
    plt.title('Mode {}'.format(Index[i]))
    plt.ylim(-31.1505, 31.1505)
    plt.xlim(50.464, 112.765)

RMS = [3.86834624,  1.86860316, 1.14843955, 1.0077887,  0.82560262]
for i in range(5):
    plt.setp(ax1[i].get_yticklabels(), visible=False)
for i in range(5):
    plt.setp(ax1[i].get_xticklabels(), visible=False)
    ax1[i].set_xlabel('{:.6f}'.format(RMS[i]))
ax1[0].set_ylabel('Clamped')

coor = np.array([])
c = 0
for j in range(3):
    for i in range(5):
        c = c + 1
        print('+--------------------------------+')
        pstring = "| Plotting mode plot    {:2d} of {:2d} |".format(
            int(c), int('15'))
        print(pstring)
        print('+--------------------------------+')
        ax2[c - 1] = plt.subplot(gs[j + 1, i])
        coor = ZModes[c]
        coor = coor.cartesianCoordinates().inUnitsOf(Ang)
        rmindex = np.array([])
        for i in range(coor.shape[0]):
            if coor[i, 2] < 10:
                rmindex = np.append(rmindex, [i])
        coor = np.delete(coor, rmindex, axis=0)
        x = coor[:, 0]
        y = coor[:, 1]
        z = (coor[:, 2]) - 11.8747
        plt.tricontourf(x, y, z, vmin=v_min, vmax=v_max, levels=levels)
        plt.ylim(-30.5601, 30.5601)
        plt.xlim(20.5441, 82.8171)

RMS = [4.95659411, 2.38986229, 1.43190479, 1.23333372,  1.0212438,
       4.35913308, 2.14893905, 1.27473199, 1.08697788,  0.90374589,
       1.97812165, 1.20437783, 0.78754608, 0.69529033,  0.63039666]
yl = ['$\epsilon = 0.31$', '$\epsilon = 1.00$', '$\epsilon = 10.0$']
rangeindex = [0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 10, 0, 0, 0, 0]
c = 0
for i in range(15):
    if i != rangeindex[i]:
        plt.setp(ax2[i].get_yticklabels(), visible=False)
    else:
        plt.setp(ax2[i].get_yticklabels(), visible=False)
        ax2[i].set_ylabel('{}'.format(yl[c]))
        c = c + 1
for i in range(15):
    plt.setp(ax2[i].get_xticklabels(), visible=False)
    ax2[i].set_xlabel('{:.6f}'.format(RMS[i]))

plt.subplots_adjust(left=0.06, bottom=None, right=0.87, top=None,
                    wspace=None, hspace=None)

cax = plt.axes([0.88, 0.1, 0.03, 0.8])
cb = plt.colorbar(ticks=[v_min, 0, v_max], cax=cax)
#cb.set_label('Out of plane displacement', fontsize=12)
cb.ax.set_yticklabels(
    ['Max\nunder\nplane', 'In\nplane', 'Max\nover\nplane'], fontsize=14)
fig.suptitle(
    'RMS values and displacement plots for different modes and substrates',
    size=18)

plt.show()

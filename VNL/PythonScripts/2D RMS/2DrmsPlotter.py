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
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit

# -------------------------------------------------------------
# Create data arrays
# -------------------------------------------------------------
myfile = np.array([])
configuration = np.array([])
dynamical_matrix = np.array([])
n_modes = np.array([])

# -------------------------------------------------------------
# Load data from dynamical matrices
# -------------------------------------------------------------
nof = 3
for i in range(nof + 1):
    if i == nof:
        break
    elif i >= 9:
        # Choose file
        myfile = np.append(myfile,
                           '{}DynamicalMatrix.hdf5'.format(i + 1))
    else:
        myfile = np.append(myfile,
                           '0{}DynamicalMatrix.hdf5'.format(i + 1))
    # Load configuration with calculator
    configuration = np.append(configuration, nlread(
        myfile[i], BulkConfiguration)[-1])

    # Load DynamicalMatrix
    dynamical_matrix = np.append(
        dynamical_matrix, nlread(myfile[i], DynamicalMatrix)[-1])

    # Vibrational state to project onto
    n_modes = np.append(n_modes, (len(configuration[i]) - len(
        dynamical_matrix[i].constraints())) * 3)

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
for i in range(nof):
    pstring = "| The File {} contains {:4d} modes |".format(
        myfile[i], int(n_modes[i]))
    print(pstring)
print('+------------------------------------------------------+')

print('+-------------------------------------------------------+')
pstring = "| The File {} contains {:4d} modes |".format(
    RFDM, int(RFn_modes))
print(pstring)
print('+-------------------------------------------------------+')

# -------------------------------------------------------------
# Create and fill data dictionaries
# -------------------------------------------------------------
qpoints = {}
frequency_list = {}
projection = {}
anti_projection = {}
RMS = {}
RFfrequency_list = {}
RFprojection = {}
RFRMS = {}

print("+=========================+")
print("|    Loading datafiles    |")
print("|-------------------------|")
with open('qpoints.pickle', 'rb') as handle:
    qpoints = pickle.load(handle)
print("| (1/5): Q-points         |")
with open('frequency_list.pickle', 'rb') as handle:
    frequency_list = pickle.load(handle)
print("| (2/5): Frequency list   |")
with open('projection.pickle', 'rb') as handle:
    projection = pickle.load(handle)
print("| (3/5): Projection       |")
with open('anti_projection.pickle', 'rb') as handle:
    anti_projection = pickle.load(handle)
print("| (4/5): Anti projection  |")
with open('RMS.pickle', 'rb') as handle:
    RMS = pickle.load(handle)
print("| (5/5): RMS              |")
print("+=========================+")

print("+=========================+")
print("|    Loading datafiles    |")
print("|-------------------------|")
with open('Referencefrequency_list.pickle', 'rb') as handle:
    RFfrequency_list = pickle.load(handle)
print("| (1/3): Frequency list   |")
with open('ReferenceProjection.pickle', 'rb') as handle:
    RFprojection = pickle.load(handle)
print("| (2/3): Projection       |")
with open('RFRMS.pickle', 'rb') as handle:
    RFRMS = pickle.load(handle)
print("| (3/3): RFRMS            |")
print("+=========================+")
# -------------------------------------------------------------
# Define colormaps
# -------------------------------------------------------------
cmap, norm = cm.get_cmap('brg'), None  # 'hot','brg','seismic'

# -------------------------------------------------------------
# Plot RMS
# -------------------------------------------------------------
ymin, ymax = -1.5, 15
fig = plt.figure(figsize=(8, 8))
gs = gridspec.GridSpec(2, 2)
gs.update(wspace=0.03)

ax1 = plt.subplot(gs[0])
coor = RFconfiguration.cartesianCoordinates().inUnitsOf(Ang)
x = coor[:, 0]
y = coor[:, 1]
tagindex = np.array(RFconfiguration.indicesFromTags(tags='FixedFinal'))
RMSscale = np.array(RFRMS[4]).flatten()

#RMSscale = RMSscale[2::3]
print(RMSscale.size)
myscale = np.ones((x.size))
j = 0
myscale[tagindex] = 0
np.set_printoptions(threshold=numpy.nan)
test = np.array([])

for i in range(myscale.size):
    if myscale[i] == 1:
        myscale[i] = RMSscale[j]
        j = j + 1
    else:
        test = np.append(test, [i])
print(test)
print(myscale)
print(RMSscale)
plt.tricontourf(x, y, myscale, 20)
plt.show()
quit()
ax2 = plt.subplot(gs[1])
myscale = {}
for i in range(nof):
    myscale[i] = projection[i]  # /numpy.max(projection)
for i in range(nof):
    # print numpy.max(projection[i])
    plt.scatter(numpy.repeat(np.array([i + 1]), n_modes[i]),
                frequency_list[i].inUnitsOf(eV).flatten(
    ) * 1000, c=myscale[i], s=15 + myscale[i] * 120, marker='o',
        edgecolor='none', cmap=cmap, norm=norm)

plt.setp(ax2.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0.07, bottom=None, right=0.88, top=None,
                    wspace=None, hspace=None)
blue_line = mlines.Line2D(
    [], [], color='blue',
    label=r'Frequency for mode $1$, $2$, $4$, $6$ and $7$')
plt.legend(handles=[blue_line], loc=8)
# colorbar
# colorbar(ticks=[-1, 0, 1], orientation='vertical')
cax = plt.axes([0.89, 0.1, 0.03, 0.8])
cb = plt.colorbar(cax=cax)
cb.set_label('projection', fontsize=12)
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

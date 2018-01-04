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
        myfile = np.append(myfile, '{}DynamicalMatrix.hdf5'.format(i + 1))
    else:
        myfile = np.append(myfile, '0{}DynamicalMatrix.hdf5'.format(i + 1))
    # Load configuration with calculator
    configuration = np.append(configuration, nlread(
        myfile[i], BulkConfiguration)[-1])

    # Load DynamicalMatrix
    dynamical_matrix = np.append(
        dynamical_matrix, nlread(myfile[i], DynamicalMatrix)[-1])

    # Vibrational state to project onto
    n_modes = np.append(n_modes, (len(configuration[i]) - len(dynamical_matrix[i].constraints())
                                  ) * 3)  # len(configuration)*3#
    # projection_vibration = numpy.zeros((n_modes),dtype=float)
# -------------------------------------------------------------
# Load data from clamped reference dynamical matrix
# -------------------------------------------------------------
RFDM = '05nmDynamicalMatrix.hdf5'
RFconfiguration = nlread(RFDM, BulkConfiguration)[-1]
RFdynamical_matrix = nlread(RFDM, DynamicalMatrix)[-1]
RFn_modes = ((len(RFconfiguration) - len(RFdynamical_matrix.constraints())
              ) * 3)

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

# quit()
# -------------------------------------------------------------
# Make projection vectors
# -------------------------------------------------------------
projection_vibration = {}
for i in range(nof):
    constrained = dynamical_matrix[i].constraints()
    tmp = numpy.zeros((n_modes[i]), dtype=float)
    tmp[2::3] = 1  # project on z-motion
    projection_vibration[i] = tmp
    # print 'Projecting on: ', projection_vibration[i], n_modes[i]

# -------------------------------------------------------------
# Set qpoint
# -------------------------------------------------------------
fractional_qpoints = [0.0, 0.0, 0.0]

# -------------------------------------------------------------
# Create data dictionaries
# -------------------------------------------------------------
qpoints = {}
frequency_list = {}
projection = {}
anti_projection = {}
RMS = {}
RFfrequency_list = {}
RFprojection = {}

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
print("| (1/2): Frequency list   |")
with open('ReferenceProjection.pickle', 'rb') as handle:
    RFprojection = pickle.load(handle)
print("| (2/2): Projection       |")
print("+=========================+")
# -------------------------------------------------------------
# Define colormaps
# -------------------------------------------------------------
cmap, norm = cm.get_cmap('brg'), None  # 'hot','brg','seismic'

# -------------------------------------------------------------
# Plot projections
# -------------------------------------------------------------
Index = np.array([0, 1, 3, 5, 6])
ymin, ymax = -0.005, 0.01
fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 6])

plt.subplot(gs[0])
myscale = RFprojection[4]  # /numpy.max(projection)
plt.scatter(numpy.repeat(np.array([1]), RFn_modes),
            RFfrequency_list[4].inUnitsOf(eV).flatten(), c=myscale,
            s=15 + myscale * 120, marker='o', edgecolor='none',
            cmap=cmap, norm=norm)
for i in range(5):
    j = Index[i]
    plt.axhline(y=RFfrequency_list[4].inUnitsOf(eV).flatten()[j])
    print(RFfrequency_list[4].inUnitsOf(eV).flatten()[j])
# Set x-ticks = the symmetry points
kticks = [1]
ticklabels = ['5nm']
xticks(kticks, ticklabels)
grid(kticks)
ylabel('$\omega$ [eV]')
#ylim(ymin, ymax)

plt.subplot(gs[1])
myscale = {}
for i in range(nof):
    myscale[i] = projection[i]  # /numpy.max(projection)
for i in range(nof):
    # print numpy.max(projection[i])
    plt.scatter(numpy.repeat(np.array([i + 1]), n_modes[i]), frequency_list[i].inUnitsOf(eV).flatten(
    ), c=myscale[i], s=15 + myscale[i] * 120, marker='o', edgecolor='none', cmap=cmap, norm=norm)
for i in range(5):
    j = Index[i]
    plt.axhline(y=RFfrequency_list[4].inUnitsOf(eV).flatten()[j])
# Set x-ticks = the symmetry points
kticks = [w * 1 for w in range(nof + 2)]
ticklabels = ['0', '0.31', '1.00', '10.0']
xticks(kticks, ticklabels)
grid(kticks)
xlabel('$\epsilon$ [8.909 meV]')
#ylim(ymin, ymax)
# colorbar
# colorbar(ticks=[-1, 0, 1], orientation='vertical')
cax = plt.axes([0.91, 0.1, 0.03, 0.8])
cb = plt.colorbar(cax=cax)
cb.ax.set_label('projection')
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

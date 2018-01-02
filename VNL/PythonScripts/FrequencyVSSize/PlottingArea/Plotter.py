# -------------------------------------------------------------
# Load libraries, ProjectedPhononBandsDisplacement and switch
# matplotlib backend
# -------------------------------------------------------------
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from NanoLanguage import *
from pylab import *
from MyAnalysisFunctions import ProjectedPhononBandsDisplacement
import pickle

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
nof = 10
for i in range(nof + 1):
    if i == nof:
        break
    elif i >= 9:
        # Choose file
        myfile = np.append(myfile, '{}nmDynamicalMatrix.hdf5'.format(i + 1))
    else:
        myfile = np.append(myfile, '0{}nmDynamicalMatrix.hdf5'.format(i + 1))
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

# Display loaded matrices and modes
print('+---------------------------------------------------------------------+')
for i in range(nof):
    pstring = "| The File {} is loaded and contains {:4d} modes |".format(
        myfile[i], int(n_modes[i]))
    print(pstring)
print('+---------------------------------------------------------------------+')
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
# Calculate projected phonon dispersion
# -------------------------------------------------------------
# Create data dictionaries
qpoints = {}
frequency_list = {}
projection = {}
anti_projection = {}
RMS = {}

for i in range(nof):
    print(
        "Calculating the projected phonon dispersion for {}".format(myfile[i]))
    T = 300. * Kelvin
    qpoints[i], frequency_list[i], projection[i], anti_projection[i], RMS[i] = ProjectedPhononBandsDisplacement(
        configuration[i], dynamical_matrix[i], fractional_qpoints, projection_vibration[i], temperature=T)
    #print(qpoints[i], frequency_list[i],projection[i], anti_projection[i], RMS[i])
print("+========================+")
print("|    Saving datafiles    |")
print("--------------------------")
with open('qpoints.pickle', 'wb') as handle:
    pickle.dump(qpoints, handle, protocol=pickle.HIGHEST_PROTOCOL)
print("| (1/5): Q-points        |")
with open('frequency_list.pickle', 'wb') as handle:
    pickle.dump(frequency_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
print("| (2/5): Frequency list  |")
with open('projection.pickle', 'wb') as handle:
    pickle.dump(projection, handle, protocol=pickle.HIGHEST_PROTOCOL)
print("| (3/5): Projection      |")
with open('anti_projection.pickle', 'wb') as handle:
    pickle.dump(anti_projection, handle, protocol=pickle.HIGHEST_PROTOCOL)
print("| (4/5): Anti projection |")
with open('RMS.pickle', 'wb') as handle:
    pickle.dump(RMS, handle, protocol=pickle.HIGHEST_PROTOCOL)
print("| (5/5): RMS             |")
print("+========================+")
quit()
# -------------------------------------------------------------
# Define colormaps
# -------------------------------------------------------------
cmap, norm = cm.get_cmap('brg'), None  # 'hot','brg','seismic'

# -------------------------------------------------------------
# Plot projections
# -------------------------------------------------------------
figure()
# plot with color and without/with variable point size
myscale = {}
for i in range(nof):
    myscale[i] = projection[i]  # /numpy.max(projection)

for i in range(nof):
    #print numpy.max(projection[i])
    plotmode = 1
    if plotmode == 0:
        scatter(numpy.repeat(np.array([i + 1]), n_modes[i]), frequency_list[i].inUnitsOf(
            eV).flatten(), c=myscale[i], s=150, marker='o', edgecolor='none', cmap=cmap, norm=norm)
    elif plotmode == 1:
        scatter(numpy.repeat(np.array([i + 1]), n_modes[i]), frequency_list[i].inUnitsOf(eV).flatten(
        ), c=myscale[i], s=15 + myscale[i] * 120, marker='o', edgecolor='none', cmap=cmap, norm=norm)
# colorbar
cb = colorbar()  # colorbar(ticks=[-1, 0, 1], orientation='vertical')
cb.set_label('projection', fontsize=20)
tick_locator = matplotlib.ticker.MaxNLocator(nbins=10, prune=None)
cb.locator = tick_locator
cb.update_ticks()
# Set x-ticks = the symmetry points
kticks = [w * 1 for w in range(12)]
ticklabels = ['%i nm' % w for w in range(11)]
xticks(kticks, ticklabels)
grid(kticks)
#ymin, ymax = 0, 0.22
# ylim(ymin,ymax)
ylabel('$\omega$ [eV]')
# plt.show()
print("+=====================================+")
print("|              Saving plots           |")
print("---------------------------------------")
print("| (1/2): FrequencyModeProjections.eps |")
savefig('FrequencyModeProjections.eps')
print("| (2/2): FrequencyModeProjections.png |")
savefig('FrequencyModeProjections.png')
print("+=====================================+")

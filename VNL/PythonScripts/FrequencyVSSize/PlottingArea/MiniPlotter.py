# -------------------------------------------------------------
# Load libraries, ProjectedPhononBandsDisplacement and switch
# matplotlib backend
# -------------------------------------------------------------
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from NanoLanguage import *
from pylab import *
from MyAnalysisFunctions import ProjectedPhononBandsDisplacement

# -------------------------------------------------------------
# Load data from dynamical matrices
# -------------------------------------------------------------
# Choose file
myfile = ('DynamicalMatrix.hdf5')
# Load configuration with calculator
configuration = nlread(myfile, BulkConfiguration)[-1]
# Load DynamicalMatrix
dynamical_matrix = nlread(myfile, DynamicalMatrix)[-1]

# Vibrational state to project onto
n_modes = (len(configuration) - len(dynamical_matrix.constraints())
                              ) * 3  # len(configuration)*3#
# projection_vibration = numpy.zeros((n_modes),dtype=float)

# Display loaded matrices and modes
print('+-----------------------------------------------------------------+')
pstring = "| The File {} is loaded and contains {:4d} modes |".format(
    myfile, int(n_modes))
print(pstring)
print('+-----------------------------------------------------------------+')

# -------------------------------------------------------------
# Make projection vectors
# -------------------------------------------------------------
constrained = dynamical_matrix.constraints()
tmp = numpy.zeros((n_modes), dtype=float)
tmp[2::3] = 1  # project on z-motion
projection_vibration = tmp
print 'Projecting on: ', projection_vibration, n_modes

# -------------------------------------------------------------
# Set qpoint
# -------------------------------------------------------------
fractional_qpoints = [0.0, 0.0, 0.0]

# -------------------------------------------------------------
# Calculate projected phonon dispersion
# -------------------------------------------------------------
print(
    "Calculating the projected phonon dispersion for {}".format(myfile))
T = 300. * Kelvin
qpoints, frequency_list, projection, anti_projection, RMS = ProjectedPhononBandsDisplacement(
    configuration, dynamical_matrix, fractional_qpoints, projection_vibration, temperature=T)
#print(qpoints, frequency_list,projection, anti_projection, RMS)

# -------------------------------------------------------------
# Define colormaps
# -------------------------------------------------------------
cmap, norm = cm.get_cmap('brg'), None  # 'hot','brg','seismic'

# -------------------------------------------------------------
# Plot projections
# -------------------------------------------------------------
figure()
# plot with color and without/with variable point size
myscale = projection  # /numpy.max(projection)
print numpy.max(projection)

plotmode = 1
if plotmode == 0:
    scatter(numpy.repeat(np.array([1]), n_modes), frequency_list.inUnitsOf(
        eV).flatten(), c=myscale, s=150, marker='o', edgecolor='none', cmap=cmap, norm=norm)
elif plotmode == 1:
    scatter(numpy.repeat(np.array([1]), n_modes), frequency_list.inUnitsOf(eV).flatten(
    ), c=myscale, s=15 + myscale * 120, marker='o', edgecolor='none', cmap=cmap, norm=norm)
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
# ymin, ymax = 0, 0.22
# ylim(ymin,ymax)
ylabel('$\omega$ [eV]')
# plt.show()
savefig('FrequencyModeProjections.eps')
savefig('FrequencyModeProjections.png')

qlist = numpy.tile(np.array([1]), (1, n_modes))
# print qlist,qlist.shape,myscale.shape,frequency_list.shape#.inUnitsOf(eV/Ang)

# -------------------------------------------------------------
# Plot RMS and frequencies
# -------------------------------------------------------------
figure()
plot(qlist, myscale, marker='o', linestyle='None')
ylabel('$\sqrt{<x^2>}$ [$\AA$]')
xlabel('Hole radia [nm]')
savefig('RMS.eps')
savefig('RMS.png')

figure()
plot(qlist, frequency_list.inUnitsOf(eV), marker='o', linestyle='None')
ylabel('$\omega$ [eV]')
xlabel('Hole radia [nm]')
savefig('Frequencies.eps')
savefig('Frequencies.png')

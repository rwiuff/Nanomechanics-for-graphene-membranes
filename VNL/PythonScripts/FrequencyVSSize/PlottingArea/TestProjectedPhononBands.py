from NanoLanguage import *
# from NL.CommonConcepts.Configurations.Utilities import cartesianToFractional,fractionalToCartesian
from pylab import *

# Choose file
# 'TestPhonons1sided_1HBN_ABstackinga=2dot49.nc'
myfile = '01nmDynamicalMatrix.hdf5'

# Load configuration with calculator
configuration = nlread(myfile, BulkConfiguration)[-1]

# Load DynamicalMatrix
dynamical_matrix = nlread(myfile, DynamicalMatrix)[-1]

# Vibrational state to project onto
n_modes = (len(configuration) - len(dynamical_matrix.constraints())
           ) * 3  # len(configuration)*3#
# projection_vibration = numpy.zeros((n_modes),dtype=float)

# Make projection vector.
constrained = dynamical_matrix.constraints()
tmp = numpy.zeros((n_modes), dtype=float)
tmp[2::3] = 1  # project on z-motion
projection_vibration = tmp
print 'Projecting on: ', projection_vibration, n_modes

# Choose q-point
myroute = ['G', 'M', 'K', 'G']
fractional_qpoints = [0.0, 0.0, 0.0]

# Calculate projected phonon dispersion
from MyAnalysisFunctions import ProjectedPhononBandsDisplacement
T = 300. * Kelvin
qpoints, frequency_list, projection, anti_projection, RMS = ProjectedPhononBandsDisplacement(
    configuration, dynamical_matrix, fractional_qpoints, projection_vibration, temperature=T)

# Calculate original dispersion
bandstructure = PhononBandstructure(
    configuration=configuration,
    route=myroute,
    points_per_segment=10,
    dynamical_matrix=dynamical_matrix,
)
bs = bandstructure
E = bs.evaluate()

# Define colormaps
cmap, norm = cm.get_cmap('brg'), None  # 'hot','brg','seismic'

# Plot bandstructure.
figure()
# plot(bs._axisGUI(), E.inUnitsOf(eV), 'k-', linewidth=2) # usual ATK dispersion
# plot(bs._axisGUI(), frequency_list.inUnitsOf(eV), 'r--', linewidth=2) # own dianonalization
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
tick_locator = matplotlib.ticker.MaxNLocator(nbins=4, prune=None)
cb.locator = tick_locator
cb.update_ticks()
# Set x-ticks = the symmetry points
kticks = [tick[0] for tick in bs._ticksGUI()]
ticklabels = [tick[1].replace('&Gamma;', '$\Gamma$')
              for tick in bs._ticksGUI()]
xticks(kticks, ticklabels)
grid(kticks)
#ymin, ymax = 0, 0.22
# ylim(ymin,ymax)
ylabel('$\omega$ [eV]')
# show()
savefig('test.png')
plt.show()
qlist = np.linalg.norm(bs.qpoints(), axis=1)
# print qlist,qlist.shape,myscale.shape,frequency_list.shape#.inUnitsOf(eV/Ang)

# Choose specific points.
i_q = 2
i_ZA = 0
i_ZOp = 5


figure()
plot(qlist, myscale, marker='o', linestyle='None')
ylabel('$\sqrt{<x^2>}$ [$\AA$]')
xlabel('q [2$\pi$/a]')
scatter(qlist[i_q], myscale[i_q, i_ZA], s=220, marker='s', edgecolor='none')
scatter(qlist[i_q], myscale[i_q, i_ZOp], s=220, marker='s', edgecolor='none')

figure()
plot(qlist, frequency_list.inUnitsOf(eV), marker='o', linestyle='None')
ylabel('$\omega$ [eV]')
xlabel('q [2$\pi$/a]')
scatter(qlist[i_q], frequency_list.inUnitsOf(eV)[
        i_q, i_ZA], s=220, marker='s', edgecolor='none')
scatter(qlist[i_q], frequency_list.inUnitsOf(eV)[
        i_q, i_ZOp], s=220, marker='s', edgecolor='none')


TempList = numpy.linspace(10.01, 1000, 21) * Kelvin
nT = len(TempList)
RMS_ZA = numpy.zeros((nT), dtype=float)
RMS_ZOp = numpy.zeros((nT), dtype=float)
for iT, T in enumerate(TempList):
    print iT, T
    qpoints, frequency_list, projection, anti_projection, RMS = ProjectedPhononBandsDisplacement(
        configuration, dynamical_matrix, fractional_qpoints, projection_vibration, temperature=T)
    RMS_ZA[iT] = RMS[i_q, i_ZA]
    RMS_ZOp[iT] = RMS[i_q, i_ZOp]

figure()
plot(TempList.inUnitsOf(Kelvin), RMS_ZA, '-or')  # marker='o',linestyle='None'
plot(TempList.inUnitsOf(Kelvin), RMS_ZOp, '-sb')  # /RMS_ZOp[0]
ylabel('$\sqrt{<x^2>}$ [$\AA$]')
xlabel('T [K]')

filesave = 'g_membrane_displacement.hdf5'
# nlsave(filesave,TempList.inUnitsOf(Kelvin))
# nlsave(filesave,RMS_ZA)
# nlsave(filesave,RMS_ZOp)


show()

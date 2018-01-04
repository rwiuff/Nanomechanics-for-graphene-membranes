# =====================================================================
#
# By Frederik Grunnet Kristensen, Christoffer Vendelbo Sorensen
# & Rasmus Wiuff, s163977@student.dtu.dk
#
# Script for creating projected mode plots and trendlines
#
# =====================================================================

# -------------------------------------------------------------
# Load libraries, ProjectedPhononBandsDisplacement and switch
# matplotlib backend for HPC compatibility
# -------------------------------------------------------------
import matplotlib.pyplot as plt
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
nof = 10
for i in range(nof + 1):
    if i == nof:
        break
    elif i >= 9:
        # Choose file
        myfile = np.append(myfile,
                           '{}nmDynamicalMatrix.hdf5'.format(i + 1))
    else:
        myfile = np.append(myfile,
                           '0{}nmDynamicalMatrix.hdf5'.format(i + 1))
    # Load configuration with calculator
    configuration = np.append(configuration, nlread(
        myfile[i], BulkConfiguration)[-1])

    # Load DynamicalMatrix
    dynamical_matrix = np.append(
        dynamical_matrix, nlread(myfile[i], DynamicalMatrix)[-1])

    # Vibrational state to project onto
    n_modes = np.append(n_modes,
                        (len(configuration[i]) -
                         len(dynamical_matrix[i].constraints())) * 3)
    # projection_vibration = numpy.zeros((n_modes),dtype=float)

# Display loaded matrices and modes
print('+-------------------------------------------------------+')
for i in range(nof):
    pstring = "| The File {} contains {:4d} modes |".format(
        myfile[i], int(n_modes[i]))
    print(pstring)
print('+-------------------------------------------------------+')

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
# Load pickled data
# -------------------------------------------------------------
# Create data dictionaries
qpoints = {}
frequency_list = {}
projection = {}
anti_projection = {}
RMS = {}
# Load data into dictionaries
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
    # print numpy.max(projection[i])
    plotmode = 1
    if plotmode == 0:
        scatter(numpy.repeat(np.array([i + 1]), n_modes[i]),
                frequency_list[i].inUnitsOf(
            eV).flatten(), c=myscale[i], s=150, marker='o',
            edgecolor='none', cmap=cmap, norm=norm)
    elif plotmode == 1:
        scatter(numpy.repeat(np.array([i + 1]), n_modes[i]),
                frequency_list[i].inUnitsOf(eV).flatten(),
                c=myscale[i], s=15 + myscale[i] * 120, marker='o',
                edgecolor='none', cmap=cmap, norm=norm)
# colorbar
cb = colorbar()  # colorbar(ticks=[-1, 0, 1], orientation='vertical')
cb.set_label('projection', fontsize=20)
tick_locator = matplotlib.ticker.MaxNLocator(nbins=10, prune=None)
cb.locator = tick_locator
cb.update_ticks()
# Set x-ticks
kticks = [w * 1 for w in range(12)]
ticklabels = ['%i nm' % w for w in range(11)]
xticks(kticks, ticklabels)
grid(kticks)

# -------------------------------------------------------------
# Fit 1/r^2 plot for mode 0
# -------------------------------------------------------------
x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y = np.zeros(10)
for i in range(nof):
    y[i] = np.sort(frequency_list[i].inUnitsOf(
        eV).flatten())[0]

# Define reciprocal function


def reciproc(x, a, n):
    return a * (1 / (x**n))


# Fit reciprocal curve
popt, pcov = curve_fit(reciproc, x, y, p0=(1, 2))
perr = np.sqrt(np.diag(pcov))
x_rec0 = np.linspace(x[0], x[-1], num=len(x) * 10)
y_rec0 = reciproc(x_rec0, *popt)

# Create legend label
rec0label = r'$a={:.3e}\pm{:.3e}$' '\n' r'$n={:.3f}\pm{:.3e}$'.format(
    popt[0], perr[0], popt[1], perr[1])

# -------------------------------------------------------------
# Fit 1/r^2 plot for mode 1
# -------------------------------------------------------------
x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y = np.zeros(10)
for i in range(nof):
    y[i] = np.sort(frequency_list[i].inUnitsOf(
        eV).flatten())[1]

# Define reciprocal function


def reciproc(x, a, n):
    return a * (1 / (x**n))


# Fit reciprocal curve
popt, pcov = curve_fit(reciproc, x, y, p0=(1, 2))
perr = np.sqrt(np.diag(pcov))
x_rec1 = np.linspace(x[0], x[-1], num=len(x) * 10)
y_rec1 = reciproc(x_rec1, *popt)

# Create legend label
rec1label = r'$a={:.3e}\pm{:.3e}$' '\n' r'$n={:.3f}\pm{:.3e}$'.format(
    popt[0], perr[0], popt[1], perr[1])

# -------------------------------------------------------------
# Show or save plots
# -------------------------------------------------------------
# Print menu
rec = r'$a \cdot \frac{1}{x^n} + b$'
showsave = 0
menu = np.array(["Show plot                  ",
                 "Show zoomed and fitted plot",
                 "Save plot                  ",
                 "Save zoomed and fitted plot"])
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
# Define graph ranges
ymin, ymax = -0.01, 0.02
# Show plot
if showsave == 1:
    ylabel('$\omega$ [eV]')
    plt.show()
# Show fitted plot
elif showsave == 2:
    ylim(ymin, ymax)
    ylabel('$\omega$ [eV]')
    ax = plt.gca()
    ax.plot(x_rec0, y_rec0, 'b-', label=rec0label)
    ax.plot(x_rec1, y_rec1, 'r-', label=rec1label)
    ax.legend(loc=3, title=rec)
    plt.show()
# Save plot
elif showsave == 3:
    ylabel('$\omega$ [eV]')
    print("+=====================================+")
    print("|              Saving plots           |")
    print("|-------------------------------------|")
    print("| (1/2): FrequencyModeProjections.eps |")
    savefig('FrequencyModeProjections.eps')
    print("| (2/2): FrequencyModeProjections.png |")
    savefig('FrequencyModeProjections.png')
    print("+=====================================+")
# Save fitted plot
elif showsave == 4:
    ylim(ymin, ymax)
    ylabel('$\omega$ [eV]')
    ax = plt.gca()
    ax.plot(x_rec0, y_rec0, 'b-', label=rec0label)
    ax.plot(x_rec1, y_rec1, 'r-', label=rec1label)
    ax.legend(loc=3, title=rec)
    print("+============================================+")
    print("|                 Saving plots               |")
    print("|--------------------------------------------|")
    print("| (1/2): FrequencyModeProjectionsZoomFit.eps |")
    savefig('FrequencyModeProjectionsZoomFit.eps')
    print("| (2/2): FrequencyModeProjectionsZoomFit.png |")
    savefig('FrequencyModeProjectionsZoomFit.png')
    print("+============================================+")

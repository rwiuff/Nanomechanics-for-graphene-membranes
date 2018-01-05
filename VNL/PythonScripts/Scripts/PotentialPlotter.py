# -------------------------------------------------------------
# Import Libraries
# -------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------
# List of Lattice Constants
# -------------------------------------------------------------
a_list = np.linspace(1, 20, num=1000)

# -------------------------------------------------------------
# Calculate total energy for each lattice constant
# -------------------------------------------------------------
Etot = np.zeros((len(a_list)), dtype=float)
for id_a, lattice_constant in enumerate(a_list):
    lattice = Hexagonal(a=lattice_constant * Angstrom, c=6.708 * Angstrom)
    elements = [Carbon] * 2
    coordinates = [(0, 0, 0), (0.33333, 0.66667, 0)]
    configuration = BulkConfiguration(lattice, elements,
                                      fractional_coordinates=coordinates)
    configuration = configuration.repeat(2, 2, 1)
    configuration = configuration.center()
    potentialSet = Tersoff_C_2010()
    calculator = TremoloXCalculator(parameters=potentialSet)
    configuration.setCalculator(calculator)
    total_energy = TotalEnergy(configuration)
    Etot[id_a] = total_energy.evaluate().inUnitsOf(eV)
plt.plot(a_list, Etot)
ymin, ymax = -70, 100
xmin, xmax = 1, 6
plt.ylim(ymin, ymax)
plt.xlim(xmin, xmax)
plt.show()

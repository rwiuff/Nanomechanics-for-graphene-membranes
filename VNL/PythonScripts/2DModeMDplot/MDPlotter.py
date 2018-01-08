# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as colors
from NanoLanguage import *
from pylab import *
import pickle
import numpy as np

# ------------------------------------------------------------
# Load configurations
# ------------------------------------------------------------
RFV = nlread('RFVconfiguration.hdf5', BulkConfiguration)[-1]
RFVeigen = nlread('RFVconfiguration.hdf5')[-1]

V = nlread('Vconfiguration.hdf5', BulkConfiguration)[-1]
Veigen = nlread('Vconfiguration.hdf5')[-1]

# ------------------------------------------------------------
# Setup simulation method
# ------------------------------------------------------------
method = Langevin(
    initial_velocity=None,
    time_step=1.0 * fs,
    reservoir_temperature=300 * Kelvin,
    friction=0 * fs**-1
)

# ------------------------------------------------------------
# Calculator for clamped system
# ------------------------------------------------------------
potentialSet = Tersoff_C_2010()
calculator = TremoloXCalculator(parameters=potentialSet)
calculator.setVerletListsDelta(0.25 * Angstrom)
RFV.setCalculator(calculator)
RFV.update()

# ------------------------------------------------------------
# Run Sumilation for clamped system
# ------------------------------------------------------------
constraints = RFV.indicesFromTags('FixedFinal')
C_md_trajectory = MolecularDynamics(
    configuration=RFV,
    constraints=constraints,
    steps=50,
    method=method,
)

ClampedCoordinates = C_md_trajectory.coordinates().inUnitsOf(Ang)
Timecoor = C_md_trajectory.times().inUnitsOf(fs)
Eigencoor = np.array([])
for i in range(ClampedCoordinates.size):
    dot = np.dot(ClampedCoordinates[i],RFVeigen[i])
    print(dot)
quit()

# ------------------------------------------------------------
# Plot simulation
# ------------------------------------------------------------

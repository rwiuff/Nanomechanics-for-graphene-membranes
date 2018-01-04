# =====================================================================
#
# By Tue Gunst (edited by Rasmus Wiuff)
#
# Script for optimising the geometry of
# a single layered membrane structure
# and define the atomic potential
#
# =====================================================================

# -------------------------------------------------------------
# Import Libraries
# -------------------------------------------------------------

from NanoLanguage import *

# -------------------------------------------------------------
# Load File
# -------------------------------------------------------------
path = 'SheetFixed.hdf5'
bulk_configuration = nlread(path, BulkConfiguration)[-1]

# -------------------------------------------------------------
# Calculator
# -------------------------------------------------------------
# Define potential for the atomic structure
potentialSet = Tersoff_C_2010()
calculator = TremoloXCalculator(parameters=potentialSet)
calculator.setVerletListsDelta(0.25 * Angstrom)

bulk_configuration.setCalculator(calculator)
nlprint(bulk_configuration)
bulk_configuration.update()

# -------------------------------------------------------------
# Optimize Geometry
# -------------------------------------------------------------
# Optimise the geometry and positions of atoms according to the
# defined potential
bulk_configuration = OptimizeGeometry(
    bulk_configuration,
    max_forces=0.001 * eV / Ang,
    max_stress=0.001 * eV / Ang**3,
    max_steps=200,
    max_step_length=0.2 * Ang,
    trajectory_filename=None,
    optimizer_method=LBFGS(),
    constrain_bravais_lattice=True,
)
nlsave('SheetRelaxed.hdf5', bulk_configuration)
nlprint(bulk_configuration)

# -------------------------------------------------------------
# Optimize Geometry
# -------------------------------------------------------------
# Optimise the geometry and positions of atoms according to the
# defined potential
# and with constrained atoms
myConstraints = bulk_configuration.indicesFromTags('FixedFinal')
bulk_configuration = OptimizeGeometry(
    bulk_configuration,
    max_forces=0.001 * eV / Ang,
    max_stress=0.001 * eV / Ang**3,
    max_steps=200,
    max_step_length=0.2 * Ang,
    constraints=myConstraints,
    trajectory_filename=None,
    optimizer_method=LBFGS(),
    constrain_bravais_lattice=True,
)
# Save configuration
nlsave('SheetRelaxed.hdf5', bulk_configuration)
nlprint(bulk_configuration)

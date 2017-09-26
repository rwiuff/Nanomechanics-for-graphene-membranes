# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
lattice = Hexagonal(3.1604*Angstrom, 25.0*Angstrom)

# Define elements
elements = [Molybdenum, Sulfur, Sulfur]

# Define coordinates
fractional_coordinates = [[ 0.333333333333,  0.666666666667,  0.5           ],
                          [ 0.666666666667,  0.333333333333,  0.4365578     ],
                          [ 0.666666666667,  0.333333333333,  0.5634422     ]]

# Set up configuration
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

# -------------------------------------------------------------
# Calculator
# -------------------------------------------------------------

potentialSet = StillingerWeber_MoS_2013()
calculator = TremoloXCalculator(parameters=potentialSet)
calculator.setVerletListsDelta(0.25*Angstrom)

bulk_configuration.setCalculator(calculator)
nlprint(bulk_configuration)
bulk_configuration.update()
nlsave('C:/Users/rwiuf/Dropbox/DTU/34029 Fagprojekt/VNL tutorial/3 - Vibrational modes/MoS2_monolayer.hdf5', bulk_configuration)

# -------------------------------------------------------------
# Optimize Geometry
# -------------------------------------------------------------
bulk_configuration = OptimizeGeometry(
    bulk_configuration,
    max_forces=0.001*eV/Ang,
    max_stress=0.001*eV/Ang**3,
    max_steps=200,
    max_step_length=0.2*Ang,
    trajectory_filename=None,
    optimizer_method=LBFGS(),
    constrain_bravais_lattice=True,
)
nlsave('C:/Users/rwiuf/Dropbox/DTU/34029 Fagprojekt/VNL tutorial/3 - Vibrational modes/MoS2_monolayer.hdf5', bulk_configuration)
nlprint(bulk_configuration)

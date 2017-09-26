# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
lattice = Hexagonal(3.09371688336*Angstrom, 25.0103649793*Angstrom)

# Define elements
elements = [Molybdenum, Sulfur, Sulfur]

# Define coordinates
fractional_coordinates = [[ 0.333333333159,  0.666666666841,  0.5           ],
                          [ 0.666666666754,  0.333333333246,  0.436383546157],
                          [ 0.666666666754,  0.333333333246,  0.563616453843]]

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
nlsave('C:/Users/rwiuf/Dropbox/DTU/34029 Fagprojekt/VNL tutorial/3 - Vibrational modes/mos2_phonons.hdf5', bulk_configuration)

# -------------------------------------------------------------
# Dynamical Matrix
# -------------------------------------------------------------
dynamical_matrix = DynamicalMatrix(
    configuration=bulk_configuration,
    repetitions=(21, 21, 1),
    atomic_displacement=0.01*Angstrom,
    acoustic_sum_rule=True,
    symmetrize=True,
    finite_difference_method=Central,
    force_tolerance=1e-08*Hartree/Bohr**2,
    processes_per_displacement=1,
    log_filename_prefix='displacement_',
    use_wigner_seitz_scheme=False,
    )
nlsave('C:/Users/rwiuf/Dropbox/DTU/34029 Fagprojekt/VNL tutorial/3 - Vibrational modes/mos2_phonons.hdf5', dynamical_matrix)

# -------------------------------------------------------------
# Phonon Density Of States
# -------------------------------------------------------------
qpoint_grid = MonkhorstPackGrid(
    na=99,
    nb=99,
    )

phonon_density_of_states = PhononDensityOfStates(
    configuration=bulk_configuration,
    dynamical_matrix=dynamical_matrix,
    qpoints=qpoint_grid,
    number_of_bands=None,
    )
nlsave('C:/Users/rwiuf/Dropbox/DTU/34029 Fagprojekt/VNL tutorial/3 - Vibrational modes/mos2_phonons.hdf5', phonon_density_of_states)
nlprint(phonon_density_of_states)

# -------------------------------------------------------------
# Phonon Bandstructure
# -------------------------------------------------------------
phonon_bandstructure = PhononBandstructure(
    configuration=bulk_configuration,
    dynamical_matrix=dynamical_matrix,
    route=['G', 'M', 'K', 'G'],
    points_per_segment=400,
    number_of_bands=All
    )
nlsave('C:/Users/rwiuf/Dropbox/DTU/34029 Fagprojekt/VNL tutorial/3 - Vibrational modes/mos2_phonons.hdf5', phonon_bandstructure)

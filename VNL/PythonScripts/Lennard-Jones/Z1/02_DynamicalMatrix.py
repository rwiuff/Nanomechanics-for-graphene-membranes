# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Load File
# -------------------------------------------------------------
path = 'SheetRelaxed.hdf5'
bulk_configuration = nlread(path, BulkConfiguration)[-1]

# -------------------------------------------------------------
# Dynamical Matrix
# -------------------------------------------------------------
myConstraints = bulk_configuration.indicesFromTags('FixedFinal')
myConstraints = map(int, myConstraints)
print(len(myConstraints))
dynamical_matrix = DynamicalMatrix(
    configuration=bulk_configuration,
    repetitions=[3, 3, 1],
    atomic_displacement=0.01 * Angstrom,
    acoustic_sum_rule=True,
    symmetrize=True,
    finite_difference_method=Central,
    force_tolerance=1e-08 * Hartree / Bohr**2,
    processes_per_displacement=1,
    log_filename_prefix='displacement_',
    use_wigner_seitz_scheme=False,
    constraints=myConstraints,
)
nlsave('01DynamicalMatrix.hdf5', dynamical_matrix)
nlsave('01DynamicalMatrix.hdf5', bulk_configuration)
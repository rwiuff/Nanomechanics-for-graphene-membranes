# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Verbosity
# -------------------------------------------------------------
verbosity = {
    'BANDSTRUCTURE': True,
    'DIAGONALIZATION_INFO': True,
    'DYNAMICAL_MATRIX': True,
    'MATRIX_SIZES': True,
    'MEMORY_INFO': True,
    'MEMORY_USAGE': True,
    'NLREAD_ERROR': True,
    'PROGRESS_BARS': True,
}
setVerbosity(verbosity)

# -------------------------------------------------------------
# Analysis from File
# -------------------------------------------------------------
path = 'DynamicalMatrix.hdf5'
configuration = nlread(path, BulkConfiguration)[-1]
dynamical_matrix = nlread(path, DynamicalMatrix)[-1]


# -------------------------------------------------------------
# Vibrational Mode
# -------------------------------------------------------------
vibrational_mode = VibrationalMode(
    configuration=configuration,
    dynamical_matrix=dynamical_matrix,
    kpoint_fractional=[0.0001, 0.0001, 0],
    mode_indices=None,
    )
nlsave('SheetVib.hdf5', vibrational_mode)

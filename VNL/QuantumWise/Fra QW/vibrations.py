# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Analysis from File
# -------------------------------------------------------------
path = 'DynamicalMatrix.hdf5'
dynamical_matrix = nlread(path, DynamicalMatrix)[-1]
configuration    = nlread(path, BulkConfiguration)[-1]

# -------------------------------------------------------------
# Vibrational Mode
# -------------------------------------------------------------
vibrational_mode = VibrationalMode(
    configuration=configuration,
    dynamical_matrix=dynamical_matrix,
    kpoint_fractional=[0, 0, 0],
    mode_indices=None,
    )
nlsave('vibrations.hdf5', vibrational_mode)

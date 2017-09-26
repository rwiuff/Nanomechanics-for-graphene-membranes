# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Analysis from File
# -------------------------------------------------------------
configuration = nlread('mos2_phonons.hdf5', object_id='BulkConfiguration_0')[0]

# -------------------------------------------------------------
# Dynamical Matrix
# -------------------------------------------------------------
dynamical_matrix = nlread('mos2_phonons.hdf5', object_id='DynamicalMatrix_0')[0]

# -------------------------------------------------------------
# Vibrational Mode
# -------------------------------------------------------------
vibrational_mode = VibrationalMode(
    configuration=configuration,
    dynamical_matrix=dynamical_matrix,
    kpoint_fractional=[0, 0, 0],
    mode_indices=[0, 1, 2, 3, 4, 5, 6, 7, 8],
    )
nlsave('mos2_vib.hdf5', vibrational_mode)

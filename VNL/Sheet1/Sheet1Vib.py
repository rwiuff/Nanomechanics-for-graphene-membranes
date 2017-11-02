# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Analysis from File
# -------------------------------------------------------------
path = 'DynamicalMatrix.hdf5'
configuration = nlread(path, object_id='BulkConfiguration_0')[0]
dynamical_matrix = nlread(path, object_id='DynamicalMatrix_0')[0]

# -------------------------------------------------------------
# Vibrational Mode
# -------------------------------------------------------------
vibrational_mode = VibrationalMode(
    configuration=configuration,
    dynamical_matrix=dynamical_matrix,
    kpoint_fractional=[0, 0, 0],
    mode_indices=None,
    )
nlsave('Sheet1Vib.hdf5', vibrational_mode)

# -------------------------------------------------------------
# Phonon Bandstructure
# -------------------------------------------------------------
phonon_bandstructure = PhononBandstructure(
    configuration=configuration,
    dynamical_matrix=dynamical_matrix,
    route=None,
    points_per_segment=20,
    number_of_bands=All
    )
nlsave('Sheet1Vib.hdf5', phonon_bandstructure)

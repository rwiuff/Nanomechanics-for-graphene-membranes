# =====================================================================
#
# By Rasmus Wiuff
#
# Script for calculating the vibrational modes for a membrane
#
# =====================================================================

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
# Save configuration
nlsave('SheetVib.hdf5', vibrational_mode)

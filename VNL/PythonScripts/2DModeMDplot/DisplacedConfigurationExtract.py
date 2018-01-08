# -------------------------------------------------------------
# Load libraries
# -------------------------------------------------------------
from NanoLanguage import *
from pylab import *
import numpy as np

# -------------------------------------------------------------
# Load data from .hdf5 files
# -------------------------------------------------------------
ZNumber = 2
if ZNumber >= 10:
    # Choose file
    myfile = '{}SheetVib.hdf5'.format(ZNumber)
else:
    myfile = '0{}SheetVib.hdf5'.format(ZNumber)

# -------------------------------------------------------------
# Load data from clamped .hdf5 file
# -------------------------------------------------------------
RFSV = '05nmSheetVib.hdf5'

# -------------------------------------------------------------
# Extract displaced configurations
# -------------------------------------------------------------
T = 300
VM = nlread(myfile, VibrationalMode)[-1]
VT = VM.movie(mode_index=2, temperature=T*Kelvin)
VE = VM.eigenvectors(mode_index=2)
VM = 0
Vconfiguration = VT.image(image_index=9)
VT = 0
nlsave('Vconfiguration.hdf5', Vconfiguration)
nlsave('Vconfiguration.hdf5', VE)

RFVM = nlread(RFSV, VibrationalMode)[-1]
RFVT = RFVM.movie(mode_index=0, temperature=T*Kelvin)
RFVE = RFVM.eigenvectors(mode_index=0)
RFVM = 0
RFVconfiguration = RFVT.image(image_index=19)
RFVT = 0
nlsave('RFVconfiguration.hdf5',RFVconfiguration)
nlsave('RFVconfiguration.hdf5',RFVE)

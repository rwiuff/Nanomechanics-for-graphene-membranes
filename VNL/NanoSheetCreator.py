# ==============================================================
# 
# Script for creating a 2 dimensional nanosheet with predefined
# hexagonal holes and sheet size using ATKPython.
#
# ==============================================================

nC = input("Repetitions along the C axis: ")
nB = input("Repetitions along the B axis: ")

sheet = NanoSheet(2,2,Carbon)

nlprint(sheet)
sheet = sheet.repeat(1,nB,nC)
sheet = sheet.center()
nlprint(sheet)

nlsave('test',sheet)
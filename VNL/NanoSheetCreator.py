# ==============================================================
#
# By Rasmus Wiuff, s163977@student.dtu.dk
# 
# Script for creating a 2 dimensional nanosheet with predefined
# hexagonal holes and sheet size using ATKPython.
#
# ==============================================================

# Import libraries
from matplotlib import path
import math

# Construct sheet
nC = input("Repetitions along the C axis: ")
nB = input("Repetitions along the B axis: ")

sheet = NanoSheet(2,2,Carbon)

sheet = sheet.repeat(1,nB,nC)
sheet = sheet.center()
sheetcoor = sheet.cartesianCoordinates()
sheetcoor = sheetcoor/Ang

str = "The sheet is {} by {} Angstrom".format(sheetcoor[-1,-1],sheetcoor[-1,-2])

print(str)
sheetcoor = numpy.array(sheetcoor)
sheetcoor = numpy.delete(sheetcoor,0,1)

# Construct Polygon
cX = input("Polygon center X-coordinate: ")
cY = input("Polygon center Y-coordinate: ")
d = input("Polygon diameter: ")
r = d/2
p1 = (cX + r, cY)
p2 = (cX + math.cos(math.pi/3)*r,cY + math.sin(math.pi/3)*r)
p3 = (cX - math.cos(math.pi/3)*r,cY + math.sin(math.pi/3)*r)
p4 = (cX - r, cY)
p5 = (cX - math.cos(math.pi/3)*r,cY - math.sin(math.pi/3)*r)
p6 = (cX + math.cos(math.pi/3)*r,cY - math.sin(math.pi/3)*r)
p = path.Path([p1,p2,p3,p4,p5,p6])

# Find atoms inside polygon
Test = p.contains_points(sheetcoor)
atoms = numpy.array([])
for i in range(Test.size):
	if Test[i] == True:
		atoms = numpy.append([atoms],[i])

# Tag the found atoms
atoms = atoms.astype(int)
atoms = atoms.tolist()


bulk_configuration = sheet

bulk_configuration.addTags('Hole', atoms)

savename = input("Input filename: ")

nlsave('savename',sheet)
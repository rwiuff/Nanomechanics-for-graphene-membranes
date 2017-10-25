# ===============================================================
#
# By Frederik Grunnet Kristensen, Christoffer Vendelbo Sorensen
# & Rasmus Wiuff, s163977@student.dtu.dk
# 
# Script for creating a 2 dimensional nanosheet with predefined
# hexagonal holes and sheet size using ATKPython.
#
# ===============================================================

# Import libraries
import numpy as np
from numpy import linalg as LA
from matplotlib import path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

# ----------------------- Construct sheet -----------------------
# Set repetitions
nA = input("Repetitions along the A axis: ")
nB = input("Repetitions along the B axis: ")

# Create lattice
lattice = Hexagonal(a=2.4612 * Angstrom,
                    c=6.709 * Angstrom)

# Set atomic elements
elements = [Carbon] * 2

# Place atoms in unit cell
coordinates = [(0, 0, 0),
               (0.33333, 0.66667, 0)]

# Create unit sheet
sheet = BulkConfiguration(lattice, elements,
                          fractional_coordinates=coordinates)

# Repeat unit sheet
sheet = sheet.repeat(nA, nB, 1)
sheet = sheet.center()

# Extract cartesian sheet coordinates
sheetcoor = sheet.cartesianCoordinates()
sheetcoor = sheetcoor / Ang

# Extract primitive vectors
pV = sheet.primitiveVectors()
pA = np.array(pV[0])
pB = np.array(pV[1])

# Show sheet size
sideA = float(LA.norm(pA) / Ang)
sideB = float(LA.norm(pB) / Ang)
size = "The graphene lattice is {:.3f} by {:.3f} Angstrom".format(
    sideA, sideB)
print(size)

# Find center coordinates
pA1 = pA / 2
pB1 = pB / 2
cX0 = pA1[0] + pB1[0]
cY0 = pA1[1] + pB1[1]

# Convert sheet coordinates to array and delete z-coordinates
sheetcoor = np.array(sheetcoor)
sheetcoor = np.delete(sheetcoor, 2, 1)

# -------------------- Construct polygon(s) ---------------------
# Loop variable and Information array
s = 1
info = np.array(["Tag name", "Center coordinate (x)",
                 "Center coordinate (y)", "Diameter", "Area"])

# Polygon creation loop
while s == 1:
    try:
        cXI = float(
            raw_input("Polygon center X-offset (default: 0): "))
    except ValueError:
        cXI = 0
    try:
        cYI = float(
            raw_input("Polygon center Y-offset (default: 0): "))
    except ValueError:
        cYI = 0
    d = float(input("Polygon diameter: "))
    r = d / 2
    cX = cX0 + cXI
    cY = cY0 + cYI
    p1 = (cX + r, cY)
    p2 = (cX + math.cos(math.pi / 3) * r,
          cY + math.sin(math.pi / 3) * r)
    p3 = (cX - math.cos(math.pi / 3) * r,
          cY + math.sin(math.pi / 3) * r)
    p4 = (cX - r, cY)
    p5 = (cX - math.cos(math.pi / 3) * r,
          cY - math.sin(math.pi / 3) * r)
    p6 = (cX + math.cos(math.pi / 3) * r,
          cY - math.sin(math.pi / 3) * r)
    p = path.Path([p1, p2, p3, p4, p5, p6, p1])

    # Find atoms inside polygon
    Test = p.contains_points(sheetcoor)
    atoms = np.array([])
    for i in range(Test.size):
        if Test[i] == True:
            atoms = np.append([atoms], [i])

    # Tag the found atoms
    atoms = atoms.astype(int)
    atoms = atoms.tolist()

    bulk_configuration = sheet

    a = (cX + math.cos(math.pi / 3) * r) - (
        cX - math.cos(math.pi / 3) * r)
    A = (3 * math.sqrt(3) * a ** 2) / 2
    area = "Area of hexagon is {:.3f} Angstrom squared".format(A)
    print(area)
    tagname = raw_input("Input tag name: ")

    bulk_configuration.addTags(tagname, atoms)

    instanceinfo = np.array(
        [tagname, "{:.3f}".format(cX), "{:.3f}".format(cY),
         "{:.3f}".format(d), "{:.3f}".format(A)])

    info = np.vstack((info, instanceinfo))

    savefig = raw_input("Save hole figure? [Y/N]: ")
    if savefig == "Y":
        fig = plt.figure()
        ax = fig.add_subplot(111)
        xp = np.transpose(sheetcoor)[0]
        yp = np.transpose(sheetcoor)[1]
        plt.plot(xp, yp, ',')
        patch = patches.PathPatch(p, facecolor='orange', lw=2)
        ax.add_patch(patch)
        ax.set_xlim(0, np.amax(sheetcoor[:, 0]))
        ax.set_ylim(0, np.amax(sheetcoor[:, 1]))
        plt.axis('equal')
        savefigfilename = raw_input(
            "Input filename for hole figure: ")
        savefigformat = "pdf"
        savefigfile = ".".join((savefigfilename, savefigformat))
        plt.savefig(savefigfile, format='pdf')
    elif savefig == "N":
        savefig = "N"

    AH = raw_input("Add hole? [Y/N]: ")
    if AH == "Y":
        s = 1
    elif AH == "N":
        s = 0

# ------------------------ Repeat sheet -------------------------
R = raw_input("Repeat sheet? [Y/N]: ")
if R == "Y":
    Areps = int(input("A-vector repetitions: "))
    sheet = sheet.repeat(Areps, 1, 1)
    Breps = int(input("B-vector repetitions: "))
    sheet = sheet.repeat(1, Breps, 1)
    # ------------------ Process repeated tags ------------------
    # Define repeater variable
    RV = int(2 * nB)
    # Define hole instances
    HI = Areps * Breps
    # Create numpy array with tag names
    tags = sheet.tags()
    tags = np.array(list(tags))
    # Create tag array
    for i in range(tags.size):
        # Tag indices for tag
        tag = sheet.indicesFromTags(tags=tags[i])
        tag = np.array(tag)
        # Sorted indices
        tag = np.sort(tag, axis=None)
        # Construct array with split tags
        tagsize = np.array(tag.shape)
        newtagsize = tagsize / HI
        newtag = np.array([])
        newtag = tag.reshape(HI, newtagsize)
        newtag = newtag.astype(int)

        for j in range(HI):
            newtagname = "".join((tags[i], str(j)))
            sheet.addTags(newtagname, newtag[j, :])

elif R == "N":
    R = "N"

savename = raw_input("Input nanosheet filename: ")
nlsave(savename, sheet)
txtformat = "txt"
savetxt = ".".join((savename, txtformat))
np.savetxt(savetxt, info, fmt='%30s %30s %30s %30s %30s',
           delimiter='|')

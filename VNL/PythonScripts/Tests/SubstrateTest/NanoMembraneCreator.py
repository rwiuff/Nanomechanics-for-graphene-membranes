# =====================================================================
#
# By Frederik Grunnet Kristensen, Christoffer Vendelbo Sorensen
# & Rasmus Wiuff, s163977@student.dtu.dk
#
# Script for creating a 2 dimensional lattice with predefined
# hexagonal holes and sheet size, using ATKPython.
#
# =====================================================================

# Import libraries
import numpy as np
from numpy import linalg as LA
from matplotlib import path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

# --------------------------- Construct sheet --------------------------
# Input user for repetitions of generating vectors
nA = input("Repetitions along the A axis: ")
nB = input("Repetitions along the B axis: ")

# Create hexagonal bravais lattice
lattice = Hexagonal(a=2.4612 * Angstrom, c=3.2 * Angstrom)

# Set atomic elements
elements = [Carbon] * 2

# Place carbon atoms in unit cell
coordinates = [(0, 0, 0), (0.33333, 0.66667, 0)]

# Create unit cell as "sheet"
sheet = BulkConfiguration(lattice, elements,
                          fractional_coordinates=coordinates)

# Repeat unit sheet along generating vectors and create substrate
sheet = sheet.repeat(nA, nB, 2)

# Increase Lattice Parameters
new_a = float(LA.norm(sheet.primitiveVectors()[0]))
new_lattice = Hexagonal(a=new_a * Angstrom, c=20 * Angstrom)
sheet.setBravaisLattice(bravais_lattice=new_lattice)
sheet = sheet.center()

# Extract cartesian sheet coordinates
sheetcoor = sheet.cartesianCoordinates()
# Remove units from coordinates
sheetcoor = sheetcoor / Ang

# Extract primitive vectors
pV = sheet.primitiveVectors()
pA = np.array(pV[0])
pB = np.array(pV[1])

# Print size of primitive vectors
sideA = float(LA.norm(pA) / Ang)
sideB = float(LA.norm(pB) / Ang)
size = "The graphene lattice is {:.3f} by {:.3f} Angstrom".format(sideA,
                                                                  sideB)
print(size)

# Find center coordinates
pA1 = pA / 2
pB1 = pB / 2
cX0 = pA1[0] + pB1[0]
cY0 = pA1[1] + pB1[1]

# Convert sheet coordinates to Numpy array and delete z-coordinates
sheetcoor = np.array(sheetcoor)
sheetcoor = np.delete(sheetcoor, 2, 1)

# ----------------------- Construct hexagon(s) -------------------------
# Loop variable and Information array
s = 1
info = np.array(
    ["Tag name", "Center coordinate (x)", "Center coordinate (y)",
     "Diameter", "Area"])

# Hexagon creation loop
while s == 1:
    try:
        # Prompt user for hexagon enter offsets
        cXI = float(
            raw_input("Hexagon center X-offset (default: 0): "))
    except ValueError:
        cXI = 0
    try:
        cYI = float(
            raw_input("Hexagon center Y-offset (default: 0): "))
    except ValueError:
        cYI = 0
    # Prompt user for hexagon diameter (greatest)
    d = float(input("Hexagon diameter: "))
    # Find radius
    r = d / 2
    # Define center coordinates
    cX = cX0 + cXI
    cY = cY0 + cYI
    # Calculate first point in hexagon
    p1 = (cX + r, cY)
    # Calculate second point in hexagon
    p2 = (cX + math.cos(math.pi / 3) * r,
          cY + math.sin(math.pi / 3) * r)
    # Calculate third point in hexagon
    p3 = (cX - math.cos(math.pi / 3) * r,
          cY + math.sin(math.pi / 3) * r)
    # Calculate fourth point in hexagon
    p4 = (cX - r, cY)
    # Calculate fifth point in hexagon
    p5 = (cX - math.cos(math.pi / 3) * r,
          cY - math.sin(math.pi / 3) * r)
    # Calculate sixth point in hexagon
    p6 = (cX + math.cos(math.pi / 3) * r,
          cY - math.sin(math.pi / 3) * r)
    # Draw path connecting the six points
    p = path.Path([p1, p2, p3, p4, p5, p6, p1])

    # Find atoms inside hexagon to Numpy array
    Test = p.contains_points(sheetcoor)
    atoms = np.array([])
    for i in range(Test.size):
        if Test[i] == True:
            atoms = np.append([atoms], [i])
            i = i + 2

    # Tag the found atoms
    atoms = atoms.astype(int)
    atoms = atoms.tolist()
    print(atoms)

    # Calculate and print area of hexagon
    a = (cX + math.cos(math.pi / 3) * r) - (
        cX - math.cos(math.pi / 3) * r)
    A = (3 * math.sqrt(3) * a ** 2) / 2
    area = "Area of hexagon is {:.3f} Angstrom squared".format(A)
    print(area)
    # Prompt user for tag name
    tagname = raw_input("Input tag name: ")
    # Add tag
    sheet.addTags(tagname, atoms)
    # Add tag information to Numpy array
    instanceinfo = np.array(
        [tagname, "{:.3f}".format(cX), "{:.3f}".format(cY),
         "{:.3f}".format(d), "{:.3f}".format(A)])

    info = np.vstack((info, instanceinfo))

    # Ask user for figure graphic
    savefig = raw_input("Save hole figure? [Y/N]: ")
    if savefig == "Y":
        # Create points with coordinates from the sheet
        fig = plt.figure()
        ax = fig.add_subplot(111)
        xp = np.transpose(sheetcoor)[0]
        yp = np.transpose(sheetcoor)[1]
        plt.plot(xp, yp, ',')
        # Create hexagonal figure
        patch = patches.PathPatch(p, facecolor='orange', lw=2)
        ax.add_patch(patch)
        ax.set_xlim(0, np.amax(sheetcoor[:, 0]))
        ax.set_ylim(0, np.amax(sheetcoor[:, 1]))
        plt.axis('equal')
        # Save figure with custom name
        savefigfilename = raw_input("Input filename for hole figure: ")
        savefigformat = "pdf"
        savefigfile = ".".join((savefigfilename, savefigformat))
        plt.savefig(savefigfile, format='pdf')
    elif savefig == "N":
        savefig = "N"
    # Prompt user to create more tags and restart loop
    AH = raw_input("Add hole? [Y/N]: ")
    if AH == "Y":
        s = 1
    elif AH == "N":
        s = 0

# Prompt user for savename
savename = raw_input("Input nanosheet filename: ")
# Save bulk configuration as .hdf5 file
nlsave(savename, sheet)
txtformat = "txt"
savetxt = ".".join((savename, txtformat))
# Create textfile with created hexagon characteristics
np.savetxt(savetxt, info, fmt='%30s %30s %30s %30s %30s',
           delimiter='|')

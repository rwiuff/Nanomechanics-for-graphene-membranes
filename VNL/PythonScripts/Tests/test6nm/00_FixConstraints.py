from NanoLanguage import *

# -------------------------------------------------------------
# Load File
# -------------------------------------------------------------
path = 'Sheet.hdf5'
bulk_configuration = nlread(path, BulkConfiguration)[-1]
constraints = bulk_configuration.indicesFromTags('Substrate')

# -------------------------------------------------------------
# Fix Constraints
# -------------------------------------------------------------

# Make indices of atoms in configuration and get constraints.
atomindices = numpy.arange(bulk_configuration.numberOfAtoms())
#constraints = constraints
# Move constrained atoms to end of list.
atomindices_new = numpy.delete(atomindices,constraints)
atomindices_new = numpy.append(atomindices_new,constraints)
constraints_new = atomindices[len(atomindices)-len(constraints)::]
new_xyz = bulk_configuration.cartesianCoordinates().inUnitsOf(Angstrom)
new_xyz = new_xyz[atomindices_new]
# Setup new configuration.
bulk_configuration._changeAtoms(positions=new_xyz*Angstrom)
constraints = constraints_new
bulk_configuration.addTags('FixedFinal',constraints)

nlsave('SheetFixed.hdf5', bulk_configuration)
nlprint(bulk_configuration)

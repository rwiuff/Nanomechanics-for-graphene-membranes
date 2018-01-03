from NanoLanguage import *

# -------------------------------------------------------------
# Load File
# -------------------------------------------------------------
path = 'SheetFixed.hdf5'
bulk_configuration = nlread(path, BulkConfiguration)[-1]

# -------------------------------------------------------------
# Lennard-Jones Calculator
# -------------------------------------------------------------
# Zeta = 0.31, 1.0, 10 for graphite, Si coupling (In SiO2) or something extremely strongly interacting, respectively.
Zeta = 0.31
# Intralayer interaction.
potentialSet_layer1 = Tersoff_C_2010(tags='Layer1')
potentialSet_layer2 = Tersoff_C_2010(tags='Layer2')
# Interlayer interaction - Lennard-Jones type.
# Define a new potential for the interlayer interaction.
lj_interlayer_potential = TremoloXPotentialSet(name="InterLayerPotential")
# Add particle type definitions for both types.
# lj_interlayer_potential.addParticleType(ParticleType.fromElement(Hydrogen))
lj_interlayer_potential.addParticleType(ParticleType.fromElement(
    Carbon, sigma=3.326 * Angstrom, epsilon=Zeta * 8.909 * meV))
# Add Lennard-Jones potentials between the carbon atoms of different layers.
lj_interlayer_potential.addPotential(
    LennardJonesPotential('C', 'C', r_cut=10.0 * Angstrom))
lj_interlayer_potential.setTags(['Layer1', 'Layer2'])

# Combine all 3 potential sets in a single calculator.
calculator = TremoloXCalculator(
    parameters=[potentialSet_layer1, potentialSet_layer2, lj_interlayer_potential])
calculator.setVerletListsDelta(0.25 * Angstrom)


bulk_configuration.setCalculator(calculator)
nlprint(bulk_configuration)
bulk_configuration.update()
myConstraints = bulk_configuration.indicesFromTags('FixedFinal')
myConstraints = map(int, myConstraints)
print(len(myConstraints))
# -------------------------------------------------------------
# Optimize Geometry
# -------------------------------------------------------------
bulk_configuration = OptimizeGeometry(
    bulk_configuration,
    max_forces=0.001 * eV / Ang,
    max_stress=0.001 * eV / Ang**3,
    max_steps=200,
    max_step_length=0.2 * Ang,
    constraints=myConstraints,
    trajectory_filename='RelaxTraj.hdf5',
    optimizer_method=LBFGS(),
    constrain_bravais_lattice=True,
)
nlsave('SheetRelaxed.hdf5', bulk_configuration)
nlprint(bulk_configuration)

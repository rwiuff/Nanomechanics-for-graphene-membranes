from NanoLanguage import *
from NL.Analysis.Mobility import fermiDistribution
#from NL.CommonConcepts.Configurations.Utilities import cartesianToFractional,fractionalToCartesian
from pylab import *
import itertools


def ProjectedPhononBandsDisplacement(configuration, dynamical_matrix, qpoints, projection_vibration, normalize_by_projection=True, temperature=300. * Kelvin):
    """ """
    """
    Private function to calculate projected phonon dispersions.
    Projects |eigenmode vectors| on input projection_vibration's (NB: displacement signs are removed).

    @param configuration            : The BulkConfiguration for which to calculate the mobility.
    @type                           : BulkConfiguration

    @param dynamical_matrix         : The dynamical matrix object.
    @type                           : DynamicalMatrix.

    @param qpoints                  : fractional q-points.

    @param projection_vibration     : The vibrational state to project onto.
    @type                           : array with length equal to the number of phonon modes.
    @default                        : No unit

    @normalize_by_projection        : Whether or not to normalize by norm of the projection vibrational state.
    @type                           : Boolean.
    @default                        : True

    # Tue Gunst, 2017
    """

    # Get dimensions
    n_modes = (len(configuration) - len(dynamical_matrix.constraints())) * 3
    n_qpoints = 1
    # Antiprojection, 1-P:
    anti_projection_vibration = numpy.ones(
        (n_modes), dtype=float) - projection_vibration

    # Calculate norm of projection vector:
    print("Calculating norm of projection vector")
    if normalize_by_projection == True:
        projection_norm = numpy.linalg.norm(projection_vibration) + 1e-12
        anti_projection_norm = numpy.linalg.norm(
            anti_projection_vibration) + 1e-12
    else:
        projection_norm = 1.0
        anti_projection_norm = 1.0

    # Masses in units of electron_mass.
    print("Calculating masses in units of electron mass")
    masses_atoms = numpy.array([e.atomicMass().inUnitsOf(
        electron_mass) for e in configuration.elements()])
    masses = numpy.repeat(masses_atoms, 3)
    masses_diagonal_matrix = numpy.diag(masses)
    # Unit factor.
    unit_factor = (1.0 * (hbar**1 * electron_mass **
                          (-0.5) * eV**(-0.5))).inUnitsOf(Ang)

    # Thermal smearing.
    thermal_smearing = (boltzmann_constant * temperature).inUnitsOf(eV)

    def coth(x):
        """
        coth(x)
        Uses numpy's tanh(x).
        """
        return 1. / tanh(x)

    # Find constraints.
    if len(dynamical_matrix.constraints()) > 0:
        constrained_displacements = numpy.sort(list(3 * dynamical_matrix.constraints()) + list(
            3 * dynamical_matrix.constraints() + 1) + list(3 * dynamical_matrix.constraints() + 2)).astype(int)
    # Calculate projections for each mode and q.
    print("Calculating projections for each mode")
    frequency_list = numpy.zeros((n_qpoints, n_modes), dtype=float)
    projection = numpy.zeros((n_qpoints, n_modes), dtype=float)
    anti_projection = numpy.zeros((n_qpoints, n_modes), dtype=float)
    RMS = numpy.zeros((n_qpoints, n_modes), dtype=float)
    qpoint = qpoints[:]


    # Calculate eigevalues and eigenvectors from dynamical matrix.
    print("Calculating eigenvalues and eigenvectors")
    eigenvalues, eigenvectors = dynamical_matrix.phononEigensystem(qpoint)
    frequency_list[0, :] = eigenvalues.inUnitsOf(eV)
    print("Done")
    for mode in numpy.arange(n_modes):
        # Remove constrained displacements.
        if len(dynamical_matrix.constraints()) > 0:
            eigenvector = numpy.delete(
                eigenvectors[:, mode], constrained_displacements)
        else:
            eigenvector = eigenvectors[:, mode]

        # Get characteristic length and mean-square mode displacement.
        eigenvalues_ph_mode = eigenvalues[mode].inUnitsOf(eV)
        inverse_characteristic_length = numpy.sqrt(2) * numpy.sqrt(abs(eigenvalues_ph_mode)) * numpy.sqrt(numpy.real(numpy.dot(
            numpy.conj(eigenvectors[:, mode].T), numpy.dot(masses_diagonal_matrix, eigenvectors[:, mode])))) / unit_factor  # 1/Ang
        ratio_thermal_frequency = abs(
            eigenvalues_ph_mode) / (2 * thermal_smearing)
        # 1./(inverse_characteristic_length**2)#
        x_squared = coth(ratio_thermal_frequency) / \
            (inverse_characteristic_length**2)

        # projection - NB: we project on absolute value of eigenmodes (direction is not important).
        CP = "({}/{}) Calculating projections for mode".format(mode+1,n_modes)
        print(CP)
        projection[0, mode] = numpy.abs(
            numpy.dot(numpy.abs(eigenvector), projection_vibration)) / projection_norm
        anti_projection[0, mode] = numpy.abs(numpy.dot(
            numpy.abs(eigenvector), anti_projection_vibration)) / anti_projection_norm
        RMS[0, mode] = numpy.sqrt(x_squared)

    # Return result
    return qpoints, frequency_list * eV, projection, anti_projection, RMS

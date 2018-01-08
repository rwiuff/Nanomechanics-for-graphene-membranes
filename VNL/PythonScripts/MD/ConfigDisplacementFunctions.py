from NanoLanguage import *
import numpy


def ChangePositionsBulk(configuration, dR):
    """
    Routine for updating configuration by dR.
    """
    R0 = configuration.cartesianCoordinates()  # .inUnitsOf(Angstrom)
    R = R0 + dR
    ConfigBulk_tmp = configuration._fastClone(cartesian_coordinates=R)
    calculator = configuration.calculator()
    ConfigBulk_tmp.setCalculator(calculator)
    return ConfigBulk_tmp


def myPhononModeDisplacedConfig(configuration, vibrational_frequencies, vibrational_states, phonon_mode, displacement_fraction=1.0, temperature=300):
    """ """
    """
    Private function to displace configuration "displacement_fraction" along a phonon mode.

    @param configuration            : The BulkConfiguration for which to calculate the mobility.
    @type                           : DeviceConfiguration

    @param temperature              : The temperature for the thermal smearing of the Bose and Fermi distributions.
    @type                           : PhysicalQuantity in units of temperature.
    @default                        : 300 Kelvin

    @vibrational_frequencies        : The vibrational frequencies.
    @vibrational_states             : The eigenmodes.
    @param phonon_mode              : The index of the mode to displace along.
    @param displacement_fraction    : The fraction of the characteristic length to displace.
    @default                        : 1.0 (can be positive or negative).

    """

    # Masses in units of electron_mass.
    masses_atoms = numpy.array([e.atomicMass().inUnitsOf(
        electron_mass) for e in configuration.elements()])
    masses = numpy.repeat(masses_atoms, 3)

    # Dimensions.
    # (len(configuration.centralRegion()) - len(dynamical_matrix.constraints()))*3
    number_of_modes = vibrational_states.shape[0]
    NA = configuration.numberOfAtoms()

    # Unit factor.
    unit_factor = (1.0 * (hbar**1 * electron_mass **
                          (-0.5) * eV**(-0.5))).inUnitsOf(Ang)

    # Thermal smearing.
    thermal_smearing = (boltzmann_constant *
                        temperature * Kelvin).inUnitsOf(eV)

    # Calculate mass-frequency scalefactor.
    masses_diagonal_matrix = numpy.diag(masses)

    # Calculate inverse characteristic length.
    mode = phonon_mode
    eigenvectors = numpy.real(vibrational_states)
    eigenvalues_ph_mode = vibrational_frequencies[mode].inUnitsOf(eV)
    inverse_characteristic_length = numpy.sqrt(2) * numpy.sqrt(abs(eigenvalues_ph_mode)) * numpy.sqrt(numpy.real(
        numpy.dot(numpy.conj(eigenvectors[:, mode].T), numpy.dot(masses_diagonal_matrix, eigenvectors[:, mode])))) / unit_factor
    # Total occupation at finite temperature and voltage.
    from NL.Analysis.Mobility import numpyBose
    bose_ph = numpyBose(eigenvalues_ph_mode, thermal_smearing)
    N = bose_ph
    sigma = numpy.sqrt(2. * N + 1) / inverse_characteristic_length
    # Get displaced configuration.
    dRt = numpy.zeros((NA * 3), dtype=float)
    # Mode and frequency of phonon mode.
    vibrational_state = numpy.real(vibrational_states[:, mode])
    #eigenvalues_ph_mode = vibrational_frequencies[mode]
    # Find index of first nonzero mode displacement.
    abs_vibrational_state = numpy.abs(vibrational_state)
    small_value = 1e-2 / number_of_modes
    index_first_nonzero_value = numpy.argmax(
        abs_vibrational_state >= small_value)  # 1e-4#1e-2/number_of_modes
    if index_first_nonzero_value == 0 and abs_vibrational_state[0] < small_value:
        print 'Error: Did not find nonzero mode displacement', myerror
    gauge_sign = numpy.sign(vibrational_state[index_first_nonzero_value])
    # update configuration.
    dRt = displacement_fraction * gauge_sign * \
        sigma * vibrational_state  # (-1)**(mode-1)* # S
    # Reshape, constrain, and employ units.
    dRt = numpy.reshape(dRt, (NA, 3))
    #dRt[iConstrainedAtoms] = 0.0
    dRt = dRt * Angstrom

    # Update configuration.
    configuration_new = ChangePositionsBulk(configuration, dRt)
    # configuration_new=ChangePositionsNoReuse(configuration,dRt)

    # Collect mode-results.
    inverse_characteristic_length = inverse_characteristic_length * 1 / Ang
    sigma = sigma * 1 / Ang

    return inverse_characteristic_length, sigma, configuration_new, dRt

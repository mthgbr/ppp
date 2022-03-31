######################################
###    AUTHOR : MATHIEU HUSSONG    ###
###   ENAC TELECOM LAB COPYRIGHT   ###
###        DATE : 16/02/2022       ###
###          VERSION : 1           ###
######################################

# THIS CODE INTERPOLATES THE POSITION AND THE TIME BIAS OF THE SATELLITES BETWEEN KNOWN POSITION AND TIME


#######################
### LIBRARY IMPORTS ###
#######################

import numpy as np
import matplotlib.pyplot as plt
from decimal import *


#################
### CONSTANTS ###
#################

C = 299792458  # speed of light in vacuum (in m/s)
DEG2RAD = np.pi / 180  # from degrees to radians conversion
ASTRO_UNIT = 149597870691  # mean distance between the Earth and the Sun (in meters)
OMEGA_EARTH = 7.2921150e-5  # rotation rate of the Earth (in rad/s)
WGS84_A = 6378137
WGS84_B = 6356752.31424518
WGS84_F = 1 / 298.257223563
WGS84_E = np.sqrt((WGS84_A**2 - WGS84_B**2) / WGS84_A**2)
WGS84_E_PRIME = np.sqrt((WGS84_A**2 - WGS84_B**2) / WGS84_B**2)
getcontext().prec = 21  # increases the numerical precision of float computations when needed

GPS = 1  # enables the pseudodistance generation with GPS satellites
GAL = 1  # enables the pseudodistance generation with Galileo satellites
GLO = 1  # enables the pseudodistance generation with Glonass satellites
BDS = 1  # enables the pseudodistance generation with Beidou satellites


###############
### CLASSES ###
###############

class SatPosition:
    """SatPosition class contains the data from the SP3 files describing
    the ECEF location of the satellites at given epochs"""
    def __init__(self):
        self.epoch_start = ""
        self.epoch_end = ""
        self.epochs = []
        self.nb_epochs = 0
        self.nb_pos = 0
        self.label = "unnamed"

    def __repr__(self):
        return "The SatPosition class '{}' contains {} satellite positions over {} epochs from the {}.".format(
            self.label, self.nb_pos, self.nb_epochs, self.epoch_start)


#################
### FUNCTIONS ###
#################

def CoM2ECEF(sun_position_ecef, satellite_position_ecef):
    """This function computes the rotation matrix from the CoM satellite frame to the ECEF frame

    %   INPUTS :
    %       Sun_position_ecef : the Sun position in the ECEF frame
    %       satellite_position_ecef : the satellite position in the ECEF frame

    %   OUTPUT :
    %       com2ecef : the rotation matrix between the CoM frame to the ECEF frame"""
    distance_earth_satellite = np.linalg.norm(satellite_position_ecef)
    sat_pos_unit_vector = [-satellite_position_ecef[0] / distance_earth_satellite,
              -satellite_position_ecef[1] / distance_earth_satellite,
              -satellite_position_ecef[2] / distance_earth_satellite]
    sun_sat_vector = sun_position_ecef - satellite_position_ecef
    e_sy_ki = np.array(np.cross(sat_pos_unit_vector, sun_sat_vector))
    norm_e_sy_ki = np.linalg.norm(e_sy_ki)
    e_sy_k = np.array(e_sy_ki / norm_e_sy_ki)
    e_sx_k = np.array(np.cross(e_sy_k, sat_pos_unit_vector))
    com2ecef = np.transpose([e_sx_k, e_sy_k, sat_pos_unit_vector])
    return np.array(com2ecef)


def ECI2ECEF(eci_position, time_univ):
    """This function converts coordinates from the ECI frame into the ECEF frame.
    A python built-in function exists and can be used instead if the
    adequate toolbox is usable.
    This function do not take into account the nutation and precession effects.

    %   INPUTS :
    %       eci_position : position in ECI to convert in ECEF (in meters)
    %       time_univ : universal time (in seconds)

    %   OUTPUT :
    %       ecef_position : position in the ECEF-frame (in meters)"""
    theta = 7.2921151467e-5 * time_univ
    ecef2eci_rotmat = [[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]]
    ecef_position = np.dot(np.linalg.inv(ecef2eci_rotmat), eci_position)
    return ecef_position


def gregorian2GPStime(year, month, day, hour, min, sec, leap_seconds=18):
    """this function returns the GPS time based on a gregorian date
    INPUTS :  a date in the gregorian calendar, given by its year, month, day, hour, minute, and second
	          leap_seconds : the number of GPS leap seconds at the given gregorian date (default = 18)

    OUTPUTS : nos = GPS time
              tod = number of seconds past last midnight"""

    hour_decimel = hour + min / 60 + sec / 3600
    if month <= 2:
        y = year - 1
        m = month + 12
    else:
        y = year
        m = month
    j = np.floor(365.25 * y)
    j = j + np.floor(30.6001 * (m + 1))
    j = j + day
    j = j + hour_decimel / 24
    j = j + 1720981.5
    num_seconds_day = 86400
    num_seconds_week = 86400 * 7
    a = np.floor(j + 0.5)
    b = a + 1537
    c = np.floor((b - 122.1) / 365.25)
    d = np.floor(365.25 * c)
    e = np.floor((b - d) / 30.6001)
    D = b
    D = D - d
    D = D - np.floor(30.6001 * e)
    D = D + (j + 0.5) % 1
    v = np.floor(j + 0.5)
    N = np.mod(v, 7)
    GPS_week = np.floor((j - 2444244.5) / 7)
    sow = (N + 1 + D%1) * num_seconds_day
    sow = np.mod(sow, num_seconds_week)
    nos = (GPS_week * num_seconds_week) + sow + leap_seconds
    tod = np.mod(nos, num_seconds_day)
    return nos, tod


def ECEF2LLA(pos_ECEF):
    """This function converts the ECEF position into the LLA one.
    %   INPUTS :
    %       position : the ECEF coordinates (in meter) of the position to convert
    %   OUTPUTS : a 1x3 vector [phi, lamda, h] with the LLA position
    %       phi : the latitude (in rad) of the position in LLA
    %       lamda : the longitude (in rad) of the position in LLA
    %       h : the altitude above the sea level (in meter) of the position"""
    p = np.sqrt(pos_ECEF[0]**2 + pos_ECEF[1]**2)
    if not p:
        lamda = 0
        phi = np.sign(pos_ECEF[2]) * np.pi / 2
        h = np.abs(pos_ECEF[2]) - WGS84_B
    else:
        theta = np.arctan(pos_ECEF[2] * WGS84_A / p / WGS84_B)
        if not pos_ECEF[0]:
            lamda = np.sign(pos_ECEF[1]) * np.pi / 2
        else:
            lamda = np.arctan(pos_ECEF[1] / pos_ECEF[0])
        phi = np.arctan((pos_ECEF[2] + WGS84_E_PRIME**2 * WGS84_B * np.sin(theta)**3) /
                        (p - WGS84_E**2 * WGS84_A * np.cos(theta)**3))
        n = WGS84_A / np.sqrt(1 - WGS84_E ** 2 * np.sin(phi) ** 2)
        if not np.cos(phi):
            h = np.abs(pos_ECEF[2]) - WGS84_B
        else:
            h = p / np.cos(phi) - n
    return np.array([phi, lamda, h])


def ECEF2ENU(beacon_position, phi_reference, lamda_reference):
    """This function converts the ECEF position into an ENU position. The
    conversion is done with respect to a reference point (needed for the ENU
    conversion).

    %   INPUTS :
    %       beacon_position : the 3x1 ECEF position to convert into the ENU frame (in meters)
    %       phi_reference : the latitude (in radians) of the reference point
    %       lambda_reference : the longitude (in radians) of the reference point

    %   OUTPUT :
    %       enu : the 3x1 position in the ENU frame (whose origin is the center of
    %           the Earth and whose axes point towards the reference point, the
    %           east of the reference point, and the North of the reference
    %           point) in meters"""
    rotation_matrix = np.zeros((3, 3))
    rotation_matrix[0, :] = [-np.sin(lamda_reference), np.cos(lamda_reference), 0]
    rotation_matrix[1, :] = [-np.sin(phi_reference)*np.cos(lamda_reference),
                             -np.sin(phi_reference)*np.sin(lamda_reference), np.cos(phi_reference)]
    rotation_matrix[2, :] = [np.cos(phi_reference) * np.cos(lamda_reference),
                             np.cos(phi_reference) * np.sin(lamda_reference), np.sin(phi_reference)]

    return np.dot(rotation_matrix, beacon_position)


def ECEF2elevation_azimuth(user_position, beacon_position):
    """This function computes the azimuth and the elevation of a point given by
    its ECEF coordinates (beacon_position) with respect to a reference location given
    by user_position. The function returns the elevation, and the azimuth of
    the beacon location as seen by the user

    %   INPUTS :
    %       user_position : the ECEF-coordinates of the reference user (in meters)
    %       beacon_position : the ECEF_coordinates of the beacon / the satellite (in meters)

    %   OUTPUTS :
    %       elevation : elevation (rad) of the satellite with respect to the receiver
    %       azimuth : azimuth (rad) of the satellite with respect to the receiver
                        azimuth is the oriented trigonometric angle (not with respect to the true north)"""
    phi, lamda, alt = ECEF2LLA(user_position)
    los_vector = beacon_position - np.transpose(user_position)
    enu_position = ECEF2ENU(los_vector, phi, lamda)
    los_range = np.linalg.norm(los_vector)
    azimuth = np.arctan2(enu_position[0], enu_position[1])
    elevation = np.arccos(np.sqrt(enu_position[0] ** 2 + enu_position[1] ** 2) / los_range) * np.sign(enu_position[2])
    return elevation, azimuth


def sun_position(gps_time, leap_seconds=18):
    """This function computes the Sun position in the ECEF frame, along with the
    distance between the Sun and the Earth.

    %   INPUTS :
    %       the date in the gregorian style at which to compute the sun position
    %       the number of leap seconds at the epoch considered (by default :18)

    %   OUTPUTS :
    %       sun_position : the position of the Sun given in the ECEF-frame"""
    time_univ = gps_time - leap_seconds
    reference_time = gregorian2GPStime(2000, 1, 1, 12, 0, 0)[0]
    time_sun = (gps_time - reference_time) / 86400 / 365.25
    eps = 23.439291 - 0.0130042 * time_sun
    sine = np.sin(eps * DEG2RAD)
    cose = np.cos(eps * DEG2RAD)
    ms = 357.5277233 + 35999.05034 * time_sun
    ls = 280.460 + 36000.770 * time_sun + 1.914666471 * np.sin(ms * DEG2RAD) + 0.019994643 * np.sin(2.0 * ms * DEG2RAD)
    rs = ASTRO_UNIT * (1.000140612 - 0.016708617 * np.cos(ms * DEG2RAD) - 0.000139589 * np.cos(2.0 * ms * DEG2RAD))
    sinl = np.sin(ls * DEG2RAD)
    cosl = np.cos(ls * DEG2RAD)
    sun_posx = rs * cosl
    sun_posy = rs * cose * sinl
    sun_posz = rs * sine * sinl
    sun_posECI = [sun_posx, sun_posy, sun_posz]
    sun_position = ECI2ECEF(sun_posECI, time_univ)
    return sun_position


def read_satellite_position(sp3_filepath, gps_antex_filepath="", gal_antex_filepath="", glo_antex_filepath="",
                            bds_antex_filepath=""):
    """this functions returns the center of mass of each satellite (GPS+GAL+GLO+BDS) sampled each 30 seconds.
    The data are drawn from a SP3 file given as input.
    The center of mass can be corrected to be the center of phase of the frequency in use, with the input of an
    ANTEX file"""
    if GPS:
        gps_antex = np.zeros((32, 3))
    if GAL:
        gal_antex = np.zeros((36, 3))
    if GLO:
        glo_antex = np.zeros((27, 3))
    if BDS:
        bds_antex = np.zeros((59, 3))

    if GPS and gps_antex_filepath:
        with open(gps_antex_filepath) as antex:
            for line in antex:
                data = line.strip("\n").split()
                sat_id = int(data[0])
                x = float(data[1]) / 1000
                y = float(data[2]) / 1000
                z = float(data[3]) / 1000
                gps_antex[sat_id-1, :] = [x, y, z]
    if GAL and gal_antex_filepath:
        with open(gal_antex_filepath) as antex:
            for line in antex:
                data = line.strip("\n").split()
                sat_id = int(data[0])
                x = float(data[1]) / 1000
                y = float(data[2]) / 1000
                z = float(data[3]) / 1000
                gal_antex[sat_id-1, :] = [x, y, z]
    if GLO and glo_antex_filepath:
        with open(glo_antex_filepath) as antex:
            for line in antex:
                data = line.strip("\n").split()
                sat_id = int(data[0])
                x = float(data[1]) / 1000
                y = float(data[2]) / 1000
                z = float(data[3]) / 1000
                glo_antex[sat_id-1, :] = [x, y, z]
    if BDS and bds_antex_filepath:
        with open(bds_antex_filepath) as antex:
            for line in antex:
                data = line.strip("\n").split()
                sat_id = int(data[0])
                x = float(data[1]) / 1000
                y = float(data[2]) / 1000
                z = float(data[3]) / 1000
                bds_antex[sat_id-1, :] = [x, y, z]

    with open(sp3_filepath) as sp3:
        satellite_positions = SatPosition()
        satellite_positions.label = sp3_filepath[:-4]
        enabled_constellations = []
        if GPS:
            enabled_constellations.append("G")
        # if GAL:
        #     enabled_constellations.append("E")
        # if GLO:
        #     enabled_constellations.append("R")
        # if BDS:
        #     enabled_constellations.append("C")
        deja_vu_epoch = False

        for line in sp3:
            data = line.strip("\n").split()

            if data[0][0] == "*":
                year = int(data[1])
                month = int(data[2])
                day = int(data[3])
                hour = int(data[4])
                minute = int(data[5])
                second = float(data[6])
                if not deja_vu_epoch:
                    satellite_positions.epoch_start = "{}/{}/{} at {}::{}::{}".format(day, month, year, hour, minute,
                                                                                    int(second))
                    deja_vu_epoch = True
                gps_time = round(gregorian2GPStime(year, month, day, hour, minute, second)[0])
                satellite_positions.epochs.append({"time": gps_time})

                satellite_positions.nb_epochs += 1

            if data[0][0] == "P" and data[0][1] in enabled_constellations:
                sat_id = data[0][1:4]
                x = float(data[1]) * 1000
                y = float(data[2]) * 1000
                z = float(data[3]) * 1000
                time = float(data[4]) * C * 1e-6
                satellite_positions.epochs[-1][sat_id] = [x, y, z, time]
                satellite_positions.nb_pos += 1

        satellite_positions.epoch_end = "{}/{}/{} at {}:{}:{}".format(day, month, year, hour, minute, int(second))

    # refinement of the satellite positions based on the ANTEX data
    for epoch in range(satellite_positions.nb_epochs):
        sun_pos_ecef = sun_position(satellite_positions.epochs[epoch]["time"])
        for sat, sat_position in satellite_positions.epochs[epoch].items():
            if sat[0] == "G":
                com2ecef = CoM2ECEF(sun_pos_ecef, sat_position[1:4])
                pco_sat = np.dot(com2ecef, gps_antex[int(sat[1:])-1, :])
                calibrated_sat_position = sat_position[1:4] + pco_sat
                satellite_positions.epochs[epoch][sat][1:4] = calibrated_sat_position
            # if sat[0] == "E":
            #     com2ecef = CoM2ECEF(sun_pos_ecef, sat_position[1:4])
            #     pco_sat = np.dot(com2ecef, gal_antex[int(sat[1:])-1, :])
            #     calibrated_sat_position = sat_position[1:4] + pco_sat
            #     satellite_positions.epochs[epoch][sat][1:4] = calibrated_sat_position
            # if sat[0] == "R":
            #     com2ecef = CoM2ECEF(sun_pos_ecef, sat_position[1:4])
            #     pco_sat = np.dot(com2ecef, glo_antex[int(sat[1:])-1, :])
            #     calibrated_sat_position = sat_position[1:4] + pco_sat
            #     satellite_positions.epochs[epoch][sat][1:4] = calibrated_sat_position
            # if sat[0] == "C":
            #     com2ecef = CoM2ECEF(sun_pos_ecef, sat_position[1:4])
            #     pco_sat = np.dot(com2ecef, bds_antex[int(sat[1:])-1, :])
            #     calibrated_sat_position = sat_position[1:4] + pco_sat
            #     satellite_positions.epochs[epoch][sat][1:4] = calibrated_sat_position

    return satellite_positions


def least_square_interpolation(satellite_positions, epoch, single_sat="", interp_degree=10, points_of_interest=20,
                               recursion_bool=True):
    """This function interpolates the satellite positions between known positions and time.
    The interpolation method consists in a weighted least square estimation of a satellite polynomial trajectory.
    INPUTS:
        satellite_positions: the SatPosition structure containing the data from sp3 files
        epoch :  the epoch at which to interpolate the satellite positions (in GPS time)
        single_sat : specify a single satellite of which to interpolate the position (as a string in the style 'G01')
                     leave single_sat="" if you want to interpolate the position of every satellites
        interp_degree : the degree of the polynomial that matches the satellite trajectory (default = 10)
        points_of_interest : the number of satellite known locations to use for interpolation (default = 20)
        recursion_bool : True if the function needs to be recursively called to improve the accuracy, False otherwise
                         please do not modify this argument, unless you know what you're doing
    OUTPUTS:
        satellite_interpolated_positions : a dictionary with the interpolated satellite positions at the input epoch
        satellite_interpolated_velocities: a dictionary with the interpolated satellite velocities at the input epoch"""
    # FINDING THE APPROPRIATE KNOWN EPOCHS TO CONSIDER TO INTERPOLATE THE TRAJECTORY AT THE EPOCH 'epoch'
    start_of_epochs_of_interest = 0

    while epoch > satellite_positions.epochs[start_of_epochs_of_interest]["time"] \
            and start_of_epochs_of_interest < satellite_positions.nb_epochs - points_of_interest/2:
        start_of_epochs_of_interest += 1
    start_of_epochs_of_interest = max(0, start_of_epochs_of_interest - int(points_of_interest / 2))

    # AT THIS MOMENT, WE KNOW WHICH SAMPLES FROM THE SP3 FLE TO CONSIDER TO INTERPOLATE THE POSITION FROM
    # INDEED, THE VARIABLE start_of_epochs_of_interest CONTAINS THE FIRST EPOCH OF THE CONSIDERED SAMPLES
    # THE N OTHER SAMPLES TO CONSIDER ARE THE N SAMPLES THAT FOLLOWS THE ONE AT EPOCH=start_of_epochs_of_interest

    satellite_interpolated_positions = {}
    for sat, sat_position in satellite_positions.epochs[start_of_epochs_of_interest].items():
        if sat != "time" and (single_sat == "" or single_sat == sat):

            # CREATION OF THE MATRICES FOR LSE
            sat_pos_x = np.zeros((points_of_interest, 1))
            sat_pos_y = np.zeros((points_of_interest, 1))
            sat_pos_z = np.zeros((points_of_interest, 1))
            sat_pos_time = np.zeros((points_of_interest, 1))
            for i in range(points_of_interest):
                sat_pos_x[i] = satellite_positions.epochs[start_of_epochs_of_interest + i][sat][0]
                sat_pos_y[i] = satellite_positions.epochs[start_of_epochs_of_interest + i][sat][1]
                sat_pos_z[i] = satellite_positions.epochs[start_of_epochs_of_interest + i][sat][2]
                sat_pos_time[i] = satellite_positions.epochs[start_of_epochs_of_interest + i][sat][3]

            coeffs = np.zeros((points_of_interest, interp_degree))
            for i in range(points_of_interest):
                for j in range(interp_degree):
                    coeffs[i, j] = (satellite_positions.epochs[start_of_epochs_of_interest + i]["time"] - epoch) ** j
            coeffs_time = np.zeros((points_of_interest, 2))
            for i in range(points_of_interest):
                for j in range(2):
                    coeffs_time[i, j] = (satellite_positions.epochs[start_of_epochs_of_interest + i]["time"] - epoch)**j

            weight = np.zeros((points_of_interest, points_of_interest))
            for i in range(points_of_interest):
                weight[i, i] = 1 / (1 + abs(satellite_positions.epochs[start_of_epochs_of_interest + i]["time"] - epoch))

            # DETERMINATION OF THE INTERPOLATION POLYNOMIAL BY COMPUTING THE LSE SOLUTION
            interp_poly_coefficients_x = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(coeffs), weight),
                                                coeffs)), np.transpose(coeffs)), weight), sat_pos_x)
            interp_poly_coefficients_y = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(coeffs), weight),
                                                coeffs)), np.transpose(coeffs)), weight), sat_pos_y)
            interp_poly_coefficients_z = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(coeffs), weight),
                                                coeffs)), np.transpose(coeffs)), weight), sat_pos_z)
            interp_poly_coefficients_time = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(coeffs_time),
                                            weight), coeffs_time)), np.transpose(coeffs_time)), weight), sat_pos_time)

            # INTERPOLATION OF THE POSITION AT THE EPOCH 'epoch'
            # AS WE DEFINED OUR TIME ORIGIN AS THE EPOCH "epoch', THE SOLUTION IS THE CONSTANT POLYNOMIAL TERM
            sat_interpolated_x = interp_poly_coefficients_x[0]
            sat_interpolated_y = interp_poly_coefficients_y[0]
            sat_interpolated_z = interp_poly_coefficients_z[0]
            sat_interpolated_time = interp_poly_coefficients_time[0]

            satellite_interpolated_positions[sat] = np.array((sat_interpolated_x[0], sat_interpolated_y[0],
                                                              sat_interpolated_z[0], sat_interpolated_time[0]))

    if recursion_bool:
        satellite_interpolated_velocities = {}
        # TO COMPUTE THE VELOCITY, WE DERIVE THE SATELLITE POSITIONS AT EPOCH = epoch - 1 second.
        # AND THEN WE USE THE APPROXIMATION THAT V = D / T
        satellite_interpolated_previous_positions = least_square_interpolation(satellite_positions, epoch-1,
                                                        single_sat, interp_degree, points_of_interest, False)
        for sat, sat_position in satellite_positions.epochs[start_of_epochs_of_interest].items():
            if sat != "time" and (single_sat == "" or single_sat == sat):
                satellite_interpolated_velocities[sat] = satellite_interpolated_positions[sat][0:3] - \
                                                         satellite_interpolated_previous_positions[sat][0:3]

        return satellite_interpolated_positions, satellite_interpolated_velocities

    return satellite_interpolated_positions


def neville_interpolation(satellite_positions, epoch, single_sat="", leap_seconds=18):
    """this function interpolates the satellite positions between known positions and time.
    The interpolation method consists in a Neville interpolation of the satellite trajectories.
    INPUTS:
        satellite_positions: the SatPosition structure containing the data from sp3 files
        epoch :  the epoch at which to interpolate the satellite positions (in GPS time)
        single_sat : specify a single satellite of which to interpolate the position (as a string in the style 'G01')
                     leave single_sat="" if you want to interpolate the position of every satellites
        leap_seconds : the number of GPS leap seconds at the time the interpolation is performed (default=18)
    OUTPUTS:
        satellite_interpolated_positions : a dictionary with the interpolated satellite positions at the input epoch
        satellite_interpolated_velocities: a dictionary with the interpolated satellite velocities at the input epoch"""
    satellite_interpolated_positions = {}
    satellite_interpolated_velocities = {}

    # FOR EVERY SATELLITE POSITION TO INTERPOLATE, WE RUN A NEVILLE INTERPOLATION METHOD
    for sat, sat_position in satellite_positions.epochs[0].items():
        if sat != "time" and (single_sat == "" or single_sat == sat):
            # tr_time represents the number of seconds since midnight
            tr_time = (epoch - leap_seconds) % 86400
            sat_position, sat_velocity, _ = neville(satellite_positions, sat, tr_time)
            satellite_interpolated_positions[sat] = sat_position
            satellite_interpolated_velocities[sat] = sat_velocity

    return satellite_interpolated_positions, satellite_interpolated_velocities


def neville(satellite_positions, id_sat, tr_time):
    """this function executes a minimization resolution based on the Neville algorithm"""

    sat_psinter_nump = 11  # interpolation degree for the position
    sat_csckinter_nump = 3  # interpolation degree for the clock offset
    Precise_time = np.linspace(0, satellite_positions.nb_epochs * 300, satellite_positions.nb_epochs + 1)
    inter_sec_delta = 300

    time_interval_ind = int(max(0, np.floor(tr_time / inter_sec_delta)))
    Xdata_comp = [float(tr_time)]

    # X estimation

    y_absp = []
    for i in range(satellite_positions.nb_epochs):
        y_absp.append(satellite_positions.epochs[i][id_sat][0])

    if time_interval_ind + sat_psinter_nump > len(Precise_time) - 1:
        x_abs = Precise_time[time_interval_ind - sat_psinter_nump:time_interval_ind + 1]
        y_abs = y_absp[time_interval_ind - sat_psinter_nump:time_interval_ind + 1]

    else:
        if not time_interval_ind:
            x_abs = Precise_time[time_interval_ind:time_interval_ind + sat_psinter_nump + 1]
            y_abs = y_absp[time_interval_ind:time_interval_ind + sat_psinter_nump + 1]
        else:
            x_abs = Precise_time[time_interval_ind - 1:time_interval_ind + sat_psinter_nump]
            y_abs = y_absp[time_interval_ind - 1:time_interval_ind + sat_psinter_nump]

    y_inter = np.zeros((len(Xdata_comp), 1))
    yveli_inter = np.zeros((len(Xdata_comp), 1))

    n = len(x_abs)

    for k in range(len(Xdata_comp)):
        xd = np.zeros(n)
        for i in range(n):
            xd[i] = abs(x_abs[i] - Xdata_comp[k])
        i = sorted(range(len(xd)), key=lambda m: xd[m])

        x_abs_sorted = []
        y_abs_sorted = []
        for ii in range(len(i)):
            x_abs_sorted.append(x_abs[i[ii]])
            y_abs_sorted.append(y_abs[i[ii]])

        P = np.zeros((n, n))
        P[:,0] = y_abs

        for i in range(n-1):
            for j in range(n-i-1):
                P[j, i + 1] = ((Xdata_comp[k] - x_abs[j]) * P[j + 1, i] + (x_abs[j+i+1]- Xdata_comp[k]) * P[j, i]) / (
                x_abs[j+i+1] - x_abs[j])

        y_inter[k] = P[0, n-1]

        D = np.zeros((n, n))
        D[:, 0] = y_abs

        for i in range(n-1):
            D[i, 1] = (D[i+1, 0] - D[i, 0]) / (x_abs[i+1] - x_abs[i])

        for i in range(1, n):
            for j in range(n-i-1):
                D[j, i + 1] = (P[j + 1, i]+ (Xdata_comp[k] - x_abs[j]) * D[j + 1, i] - P[j, i] + (
                x_abs[j + i + 1] - Xdata_comp[k]) * D[j, i]) / (x_abs[j + i + 1] - x_abs[j])

        yveli_inter[k] = D[0, n-1]

    X_inter = y_inter[0][0]
    VX_inter = yveli_inter[0][0]

    # Y estimation

    y_absp = []
    for i in range(satellite_positions.nb_epochs):
        y_absp.append(satellite_positions.epochs[i][id_sat][1])

    if time_interval_ind + sat_psinter_nump > len(Precise_time) - 1:
        x_abs = Precise_time[time_interval_ind - sat_psinter_nump + 1:time_interval_ind + 2]
        y_abs = y_absp[time_interval_ind - sat_psinter_nump + 1:time_interval_ind + 2]

    else:
        if not time_interval_ind:
            x_abs = Precise_time[time_interval_ind + 1:time_interval_ind + sat_psinter_nump + 2]
            y_abs = y_absp[time_interval_ind + 1:time_interval_ind + sat_psinter_nump + 2]
        else:
            x_abs = Precise_time[time_interval_ind:time_interval_ind + sat_psinter_nump + 1]
            y_abs = y_absp[time_interval_ind:time_interval_ind + sat_psinter_nump + 1]

    y_inter = np.zeros((len(Xdata_comp), 1))
    yveli_inter = np.zeros((len(Xdata_comp), 1))

    n = len(x_abs)

    for k in range(len(Xdata_comp)):
        xd = np.zeros(n)
        for i in range(n):
            xd[i] = abs(x_abs[i] - Xdata_comp[k])
        i = sorted(range(len(xd)), key=lambda m: xd[m])

        x_abs_sorted = []
        y_abs_sorted = []
        for ii in range(len(i)):
            x_abs_sorted.append(x_abs[i[ii]])
            y_abs_sorted.append(y_abs[i[ii]])

        P = np.zeros((n, n))
        P[:, 0] = y_abs

        for i in range(n - 1):
            for j in range(n - i - 1):
                P[j, i + 1] = ((Xdata_comp[k] - x_abs[j]) * P[j + 1, i] + (x_abs[j + i + 1] - Xdata_comp[k]) * P[
                    j, i]) / (
                                      x_abs[j + i + 1] - x_abs[j])

        y_inter[k] = P[0, n - 1]

        D = np.zeros((n, n))
        D[:, 0] = y_abs

        for i in range(n - 1):
            D[i, 1] = (D[i + 1, 0] - D[i, 0]) / (x_abs[i + 1] - x_abs[i])

        for i in range(1, n - 1):
            for j in range(n - i - 1):
                D[j, i + 1] = (P[j + 1, i] + (Xdata_comp[k] - x_abs[j]) * D[j + 1, i] - P[j, i] + (
                        x_abs[j + i + 1] - Xdata_comp[k]) * D[j, i]) / (x_abs[j + i + 1] - x_abs[j])

        yveli_inter[k] = D[0, n - 1]

    Y_inter = y_inter[0][0]
    VY_inter = yveli_inter[0][0]

    # Z estimation

    y_absp = []
    for i in range(satellite_positions.nb_epochs):
        y_absp.append(satellite_positions.epochs[i][id_sat][2])

    if time_interval_ind + sat_psinter_nump > len(Precise_time) - 1:
        x_abs = Precise_time[time_interval_ind - sat_psinter_nump + 1:time_interval_ind + 2]
        y_abs = y_absp[time_interval_ind - sat_psinter_nump + 1:time_interval_ind + 2]

    else:
        if not time_interval_ind:
            x_abs = Precise_time[time_interval_ind + 1:time_interval_ind + sat_psinter_nump + 2]
            y_abs = y_absp[time_interval_ind + 1:time_interval_ind + sat_psinter_nump + 2]
        else:
            x_abs = Precise_time[time_interval_ind:time_interval_ind + sat_psinter_nump + 1]
            y_abs = y_absp[time_interval_ind:time_interval_ind + sat_psinter_nump + 1]


    y_inter = np.zeros((len(Xdata_comp), 1))
    yveli_inter = np.zeros((len(Xdata_comp), 1))

    n = len(x_abs)

    for k in range(len(Xdata_comp)):
        xd = np.zeros(n)
        for i in range(n):
            xd[i] = abs(x_abs[i] - Xdata_comp[k])
        i = sorted(range(len(xd)), key=lambda m: xd[m])

        x_abs_sorted = []
        y_abs_sorted = []
        for ii in range(len(i)):
            x_abs_sorted.append(x_abs[i[ii]])
            y_abs_sorted.append(y_abs[i[ii]])

        P = np.zeros((n, n))
        P[:, 0] = y_abs

        for i in range(n - 1):
            for j in range(n - i - 1):
                P[j, i + 1] = ((Xdata_comp[k] - x_abs[j]) * P[j + 1, i] + (x_abs[j + i + 1] - Xdata_comp[k]) * P[
                    j, i]) / (
                                      x_abs[j + i + 1] - x_abs[j])

        y_inter[k] = P[0, n - 1]

        D = np.zeros((n, n))
        D[:, 0] = y_abs

        for i in range(n - 1):
            D[i, 1] = (D[i + 1, 0] - D[i, 0]) / (x_abs[i + 1] - x_abs[i])

        for i in range(1, n - 1):
            for j in range(n - i - 1):
                D[j, i + 1] = (P[j + 1, i] + (Xdata_comp[k] - x_abs[j]) * D[j + 1, i] - P[j, i] + (
                        x_abs[j + i + 1] - Xdata_comp[k]) * D[j, i]) / (x_abs[j + i + 1] - x_abs[j])

        yveli_inter[k] = D[0, n - 1]

    Z_inter = y_inter[0][0]
    VZ_inter = yveli_inter[0][0]

    # clock estimation

    y_absp = []
    for i in range(satellite_positions.nb_epochs):
        y_absp.append(satellite_positions.epochs[i][id_sat][3])

    if time_interval_ind + sat_csckinter_nump > len(Precise_time) - 1:
        x_abs = Precise_time[time_interval_ind - sat_csckinter_nump + 1:time_interval_ind + 2]
        y_abs = y_absp[time_interval_ind - sat_csckinter_nump + 1:time_interval_ind + 2]

    else:
        if not time_interval_ind:
            x_abs = Precise_time[time_interval_ind + 1:time_interval_ind + sat_csckinter_nump + 2]
            y_abs = y_absp[time_interval_ind + 1:time_interval_ind + sat_csckinter_nump + 2]
        else:
            x_abs = Precise_time[time_interval_ind:time_interval_ind + sat_csckinter_nump + 1]
            y_abs = y_absp[time_interval_ind:time_interval_ind + sat_csckinter_nump + 1]

    y_inter = np.zeros((len(Xdata_comp), 1))
    yveli_inter = np.zeros((len(Xdata_comp), 1))

    n = len(x_abs)

    for k in range(len(Xdata_comp)):
        xd = np.zeros(n)
        for i in range(n):
            xd[i] = abs(x_abs[i] - Xdata_comp[k])
        i = sorted(range(len(xd)), key=lambda m: xd[m])

        x_abs_sorted = []
        y_abs_sorted = []
        for ii in range(len(i)):
            x_abs_sorted.append(x_abs[i[ii]])
            y_abs_sorted.append(y_abs[i[ii]])

        P = np.zeros((n, n))
        P[:, 0] = y_abs

        for i in range(n - 1):
            for j in range(n - i - 1):
                P[j, i + 1] = ((Xdata_comp[k] - x_abs[j]) * P[j + 1, i] + (x_abs[j + i + 1] - Xdata_comp[k]) * P[
                    j, i]) / (
                                      x_abs[j + i + 1] - x_abs[j])

        y_inter[k] = P[0, n - 1]

        D = np.zeros((n, n))
        D[:, 0] = y_abs

        for i in range(n - 1):
            D[i, 1] = (D[i + 1, 0] - D[i, 0]) / (x_abs[i + 1] - x_abs[i])

        for i in range(1, n - 1):
            for j in range(n - i - 1):
                D[j, i + 1] = (P[j + 1, i] + (Xdata_comp[k] - x_abs[j]) * D[j + 1, i] - P[j, i] + (
                        x_abs[j + i + 1] - Xdata_comp[k]) * D[j, i]) / (x_abs[j + i + 1] - x_abs[j])

        yveli_inter[k] = D[0, n - 1]

    clock_inter = y_inter[0][0]
    clock_drift_inter = yveli_inter[0][0]

    Rsat = [X_inter, Y_inter, Z_inter]
    Vsat = [VX_inter, VY_inter, VZ_inter]

    relativity = (Rsat * np.transpose(Vsat))/ C**2
    clock_error = clock_inter -2 * C * relativity
    satclock_vect = [clock_error,  clock_drift_inter]
    Rsat.append(clock_inter)

    return Rsat, Vsat, satclock_vect


def lagrande_interpolation(t, v, T):
    """% This file intepolates a function defined by a set of abscissa and a
    % corresponding set of values. The interpolation is done at given points
    % thanks to the Lagrande method and the interpolated value is returned.

    %   INPUTS :
    %       t : a vector containing the abscissa points at which the function
    %           is sampled ; t = [x1 x2 ... x_n]
    %       v : a vector containing the value of the function at the points
    %           given by t ; v = [f(x1) f(x2) ... f(x_n)]
    %       T : a vector or a value of the point(s) where to interpolate the
    %           value."""
    n = len(t)
    sum = 0
    for i in range(n):
        prod = v[i]
        for j in range(n):
            if i != j:
                prod = prod * (T - t[j]) / (t[i] - t[j])
        sum += prod
    return sum


def check_interpolation_accuracy(accurate_pos_filepath, interpolated_data):
    """This function was designed to verify the interpolation consistency. It is now of no use."""
    nb_samples = len(interpolated_data[:, 0])
    accurate_data = np.zeros((1800, 3))
    time = np.linspace(1, nb_samples, nb_samples)
    with open(accurate_pos_filepath) as sat_accurate_pos:
        nb_line = 0
        for line in sat_accurate_pos:
            data = line.strip("\n").split(',')
            accurate_data[nb_line, :] = [float(data[0]), float(data[1]), float(data[2])]
            nb_line += 1
    print(accurate_data[0, :])
    plt.plot(time[:nb_samples], accurate_data[:nb_samples, 0], "-k", linewidth=2, label="ref")
    plt.plot(time[:nb_samples], interpolated_data[:, 0], "-r", linewidth=2, label="interp")
    plt.title("x coordinate over time")
    plt.xlabel("time (s)")
    plt.ylabel("x (m)")
    plt.legend()
    plt.show()
    plt.plot(time[:nb_samples], accurate_data[:nb_samples, 0] - interpolated_data[:, 0], "-g", linewidth=2)
    plt.title("x error over time")
    plt.xlabel("time (s)")
    plt.ylabel("x (m)")
    plt.show()

def euclidian_distance(user_position, satellite_positions, epoch, elevation_mask=0):
#def euclidian_distance(user_position, satellite_positions, epoch):
    """returns the Euclidian distances between the user and all the visible satellites"""

    satellite_positions_at_reception_epoch, satellite_velocities_at_reception_epoch = \
        neville_interpolation(satellite_positions, epoch)

    satellite_elevations_azimuths = {}
    for sat, sat_position in satellite_positions_at_reception_epoch.items():
        elevation, azimuth = ECEF2elevation_azimuth(user_position, sat_position[0:3])
        satellite_elevations_azimuths[sat] = [elevation, azimuth]

    euclidian_distances = {}
    for sat, sat_position in satellite_positions_at_reception_epoch.items():
        if satellite_elevations_azimuths[sat][0] * 180 / np.pi > elevation_mask:
            distance = np.sqrt((user_position[0] - sat_position[0])**2 + (user_position[1] - sat_position[1])**2
                               + (user_position[2] - sat_position[2])**2)
            euclidian_distances[sat] = {"distance": distance, "tx_time": epoch}



    # refinement of the distance and the emission time taking into account the displacement of the satellite in the
    # transmission duration :
        tolerable_error = 1/2**26 # delta_tx error = 1.5*10^-8 s ==> delta_pos error < 5*10^-5 m
    for sat, sat_dist_and_time in euclidian_distances.items():
        error = Decimal(epoch) - Decimal(sat_dist_and_time["tx_time"]) - Decimal(sat_dist_and_time["distance"]) / C
        tx_epoch = Decimal(sat_dist_and_time["tx_time"])
        sat_distance = sat_dist_and_time["distance"]
        while abs(error) > tolerable_error:
            tx_epoch += error
            pos_dict, vel_dict = neville_interpolation(satellite_positions,
                                tx_epoch, sat)
            new_sat_pos = pos_dict[sat]
            new_sat_vel = vel_dict[sat]
            sat_distance = np.sqrt((user_position[0] - new_sat_pos[0])**2 + (user_position[1] - new_sat_pos[1])**2
                               + (user_position[2] - new_sat_pos[2])**2)
            error = Decimal(epoch) - tx_epoch - Decimal(sat_distance) / C

        euclidian_distances[sat] = {"distance": sat_distance,
                            "tx_time": tx_epoch,
                            "pos": new_sat_pos[0:3], "vel": new_sat_vel, "clock": new_sat_pos[3]}

        # compensating the earth rotation at the actual emission time
        theta = OMEGA_EARTH * float((Decimal(epoch) - tx_epoch))
        rotation_theta = np.array([[np.cos(theta), np.sin(theta), 0],
                                    [-np.sin(theta), np.cos(theta), 0],
                                   [0, 0, 1]])
        corrected_sat_pos = np.transpose(np.dot(rotation_theta, [[new_sat_pos[0]], [new_sat_pos[1]], [new_sat_pos[2]]]))[0]
        corrected_sat_vel = np.transpose(np.dot(rotation_theta, [[new_sat_vel[0]], [new_sat_vel[1]], [new_sat_vel[2]]]))[0]

        corrected_distance = np.sqrt((user_position[0] - corrected_sat_pos[0])**2 + (user_position[1] - corrected_sat_pos[1])**2
                               + (user_position[2] - corrected_sat_pos[2])**2)
        euclidian_distances[sat] = {"distance": corrected_distance,
                                    "tx_time": tx_epoch,
                                    "pos": corrected_sat_pos, "vel": corrected_sat_vel, "clock": new_sat_pos[3]}

    return euclidian_distances

def run():
    # TESTS IF THE CODE IS WORKING BY GETTING HOLD OF THE SATELLITE POSITIONS AT EVERY SP3 EPOCHS
    satellite_positions = read_satellite_position("C:/Users/Levina Shewina/PycharmProjects/SATELLITE_POSITION_FILES/WUM0MGXULA_20201010000_01D_05M_ORB.SP3",
                                                  "C:/Users/Levina Shewina/PycharmProjects/ANTEX_FILES/gpsantex2020.txt",
                                                  "C:/Users/Levina Shewina/PycharmProjects/ANTEX_FILES/galileoantex2020.txt",
                                                  "C:/Users/Levina Shewina/PycharmProjects/ANTEX_FILES/glonassantex2020.txt",
                                                  "C:/Users/Levina Shewina/PycharmProjects/ANTEX_FILES/beidouantex2020.txt")
    # TESTS IF THE CODE IS WORKING BY INTERPOLATING THE POSITION AT A RANDOM EPOCH
    time = [2020, 4, 9, 1, 58, 00]
    #response = {}
    i = 0
    while i < 86400:
        epoch_to_test_the_interpolation = gregorian2GPStime(*time)[0]
        #print(epoch_to_test_the_interpolation)

        lse_positions = least_square_interpolation(satellite_positions, epoch_to_test_the_interpolation)
        neville_positions = neville_interpolation(satellite_positions, epoch_to_test_the_interpolation)
        euclidian_distances = euclidian_distance([WGS84_A, 0, 0], satellite_positions, epoch_to_test_the_interpolation)
        print("\n".join(str(time[3] * 3600 + time[4] * 60 + time[5]) + ": " + "{!r}: {!r},".format(k, v) for k, v in euclidian_distances.items()))

        time[5] = (time[5] + 1) % 60
        if time[5] == 0:
            time[4] = (time[4] + 1) % 60
            if time[4] == 0:
                time[3] = (time[3] + 1) % 3600
        i += 1

        #response[str(time[4] * 60 + time[5])] = euclidian_distances
    #return response



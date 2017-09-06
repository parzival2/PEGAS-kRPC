import numpy as np
import init_simulation
import unified_powered_flight_guidance as upfg

mu = init_simulation.mu
R = init_simulation.R
period = init_simulation.period


class LaunchSite:
    def __init__(self, latitude, longitude, altitude):
        self.latitude = float(latitude)
        self.longitude = float(longitude)
        self.altitude = float(altitude)


def launch_targeting(site, altitude, opposing_apsis, inclination, deltaLAN):
    """
    Creates a UPFG-compatible target structure as well as a launch azimuth for
    first stage guidance from a given launch site and definition of a target
    orbit.
    Assumes the insertion will occur in one apsis of the target orbit. Handles
    LAN specially - since the Earth does not rotate in the simulation, (it is
    assumed that at launch (beginning of any simulation) the 0 degree meridian
    points in a reference direction - ie. a zero-LAN orbit's ascending node
    passes right over it at this time) arbitrarily choosing the target LAN
    makes no sense. Instead, an assumption is made that the launch occurs some
    time before the launch site rotates under the target orbit. This time can
    be specified by the deltaLAN param. This function will return the LAN of a
    targeted orbit, so in order to simulate a launch into a plane with any
    desired LAN, a simple adjustment of time to launch is possible: one simply
    needs to subtract the returned LAN from their desired LAN and convert this
    to seconds (ie. calculate how long will it take until the Earth rotates by
    this amount).

    :param site: Struct defining the launch site.
    :param altitude: Desired cutoff altitude (km above sea level). Becomes one
        apsis of the target orbit.
    :param opposing_apsis: Opposing apsis of the target orbit (km above sea level).
    :param inclination: inclination of the target orbit in degrees.
    :param deltaLAN: Lift-off will occur when the launch site is that many
        degrees before rotating right under the target orbit.
        Zero means launch right when this happens. Negative values
        mean launch is late, after the target orbit passed over the
        launch site. Values of 1.5-2.5 work pretty good, depending
        on vehicle. Too little or too much can break UPFG.

    :return: (lan, azimuth, target)
        lan      LAN of the orbit passing (degrees).
        azimuth  Launch azimuth corrected for the velocity bonus from
            Earth's rotation (degrees from earth to north, CCW).
        target   Target struct, compatible with UPFG.

    Calculation of the launch azimuth. Math largely based on:
    http://www.orbiterwiki.org/wiki/Launch_Azimuth
    """
    global mu      # Global variable, standard gravity parameter of the body;
                   # gravity constant * mass of the body (kg).
    global R       # Global variable, radius of the body (m).
    global period  # Global variable, period of Earth's rotation (seconds).

    if inclination < site.latitude:
        print('Target inclination below launch site latitude. Assuming you know what youre doing - returning 0 azimuth.')
        azimuth = 0
    else:
        binertial = np.rad2deg(np.arcsin(np.cos(np.deg2rad(inclination))/np.cos(np.deg2rad(site.latitude))))  # launch azimuth with no regard for Earth rotation
        vorbit = np.sqrt(mu/(R+altitude*1000))                            # orbital velocity magnitude
        vEarthrot = (2*np.pi*R/period)*np.cos(np.deg2rad(site.latitude))  # v gained from Earth rotation
        vrotx = vorbit*np.sin(np.deg2rad(binertial))-vEarthrot
        vroty = vorbit*np.cos(np.deg2rad(binertial))
        azimuth = np.rad2deg(np.arctan2(vroty, vrotx))                    # corrected launch azimuth

    # LAN of the orbit passing directly over the launch site at lift-off.
    # Math derived from Napier's rules for spherical triangles.
    lan = site.longitude - np.rad2deg(np.arcsin(min(1, np.tan(np.deg2rad(90-inclination)) * np.tan(np.deg2rad(site.latitude)))))
    lan = np.mod(lan + 360 + deltaLAN, 360)                           # add slip and reduce result into 0-360 range

    # Target orbit normal vector obtained by rotation of an unit vector
    # using matrices of rotations (Ry is not used). Vector points south
    # instead of north (ie. [0,0,-1]) because UPFG requires it to point
    # opposite to a direction of vector of angular momentum for the orbit.
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(np.deg2rad(inclination)), -np.sin(np.deg2rad(inclination))],
                   [0, np.sin(np.deg2rad(inclination)), np.cos(np.deg2rad(inclination))]])  # about x for inclination (preserve zero node)
    Ry = np.array([[np.cos(0), 0, np.sin(0)],
                   [0, 1, 0],
                   [-np.sin(0), 0, np.cos(0)]])                                             # in case we needed it for something
    Rz = np.array([[np.cos(np.deg2rad(lan)), -np.sin(np.deg2rad(lan)), 0],
                   [np.sin(np.deg2rad(lan)), np.cos(np.deg2rad(lan)), 0],
                   [0, 0, 1]])                                                              # about z for node
    target_plane_normal = np.matmul(Rz, np.matmul(Rx, np.array([0, 0, -1])))
    # Target velocity from the vis-viva equation, step by step for clarity.
    target_velocity = 2.0 / (R+altitude*1000)
    target_velocity = target_velocity - 1.0 / (R+(altitude+opposing_apsis)*1000/2)
    target_velocity = np.sqrt(mu*target_velocity)
    # Flight path angle is always zero (ie. cutoff is desired to occur in
    # an apsis).
    target = upfg.Target(0, target_plane_normal, R+altitude*1000, target_velocity)
    return lan, azimuth, target

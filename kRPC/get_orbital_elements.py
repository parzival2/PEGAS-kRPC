import numpy as np
import init_simulation

mu = init_simulation.mu
R = init_simulation.R

# function [ap, pe, sma, ecc, inc, lan, aop, tan]
def get_orbital_elements(r, v):
    """
    Calculates complete set of Keplerian elements of a 3D orbit given by
    position and velocity vectors. Equations taken from:
    https://en.wikibooks.org/wiki/Astrodynamics/Classical_Orbit_Elements
    http://space.stackexchange.com/a/1919

    :param r: Position XYZ vector relative to the center of reference
              frame (m).
    :param v: Velocity XYZ vector (m/s).
    :return: (ap. pe, sma, ecc, inc, lan, aop, tan)
        ap         Apoapsis from the body's surface (km).
        pe         Periapsis from the body's surface (km).
        sma        Semi-major axis (m).
        ecc        Eccentricity
        inc        Inclination (deg).
        lan        Longitude of ascending node (deg).
        aop        Argument of periapsis (deg).
        tan        True anomaly (deg).
    """
    global mu   # Global variable, standard gravity parameter of the body;
                # gravity constant * mass of the body (kg).
    global R    # Global variable, radius of the body (m).
    # angular momentum
    h = np.cross(r, v)
    # ascending node
    n = np.cross(np.array([0, 0, 1]), h)
    # specific mechanical energy
    E = np.linalg.norm(v)**2/2-mu/np.linalg.norm(r)
    # semi-major axis
    sma = -0.5*mu/E
    # eccentricity vector
    e = (np.linalg.norm(v)**2/mu - 1/np.linalg.norm(r))*r - np.vdot(r, v)/mu*v
    ecc = np.linalg.norm(e)
    # inclination
    inc = np.rad2deg(np.arccos(np.vdot(np.array([0, 0, 1]), h)/np.linalg.norm(h)))
    # longitude of the ascending node
    lan = np.rad2deg(np.arccos(np.vdot(np.array([1, 0, 0]), n)/np.linalg.norm(n)))
    if n[1] < 0:
        lan = 360-lan
    # argument of periapsis
    aop = np.rad2deg(np.arccos(np.vdot(n, e)/(np.linalg.norm(n)*np.linalg.norm(e))))
    if e[2] < 0:
        aop = 360-aop
    # true anomaly
    tan = np.rad2deg(np.arccos(np.vdot(e, r)/(np.linalg.norm(e)*np.linalg.norm(r))))
    if np.vdot(r,v) < 0:
        tan = 360-tan
    # Ap, Pe
    ap = (sma*(1+ecc)-R)/1000
    pe = (sma*(1-ecc)-R)/1000
    return ap, pe, sma, ecc, inc, lan, aop, tan

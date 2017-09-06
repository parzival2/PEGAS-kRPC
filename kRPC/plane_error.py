import numpy as np


def plane_error(results, target):
    """
    Computes angle between target orbital plane and actually achieved plane.
    
    :param results: Results struct as output by flight_manager (NOT flight_sim_3d).
    :param target: Target struct as output by launch_targeting.
    :return: Angle between the two orbital planes.
    """
    inc = results.powered[results.n-1].orbit.inc
    lan = results.powered[results.n-1].orbit.lan
    
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(np.deg2rad(inc)), -np.sin(np.deg2rad(inc))],
                   [0, np.sin(np.deg2rad(inc)), np.cos(np.deg2rad(inc))]])
    Rz = np.array([[np.cos(np.deg2rad(lan)), -np.sin(np.deg2rad(lan)), 0],
                   [np.sin(np.deg2rad(lan)), np.cos(np.deg2rad(lan)), 0],
                   [0, 0, 1]])
    reached = np.matmul(Rz, np.matmul(Rx, np.array([0, 0, -1])))
    error = np.rad2deg(np.arccos(np.vdot(target.normal, reached)))
    return error

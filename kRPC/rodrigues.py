import numpy as np
from unit import unit


def rodrigues(vector, axis, angle):
    """
    Rotates a given vector about a given axis by a given angle using Rodrigues
    rotation formula:
    https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula

    :param vector: XYZ vector to be rotated.
    :param axis: XYZ vector representing the axis of rotation (will be
        normalized so its magnitude does not matter).
    :param angle: Angle of rotation (degrees).
    :return: Rotated XYZ vector.
    """
    axis = unit(axis);
    rotated = vector*np.cos(np.deg2rad(angle))
    rotated = rotated + np.cross(axis, vector)*np.sin(np.deg2rad(angle))
    rotated = rotated + axis*np.vdot(axis,vector)*(1-np.cos(np.deg2rad(angle)))
    return rotated

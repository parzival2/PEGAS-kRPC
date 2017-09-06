import numpy as np
from unit import unit


def get_angle_from_frame(vector, frame, type):
    """
    Calculates pitch or yaw angle of a vector in a given reference frame.
       
    :param vector: XYZ vector 
    :param frame: 3x3 matrix of unit basis vectors given row-wise in the
        following order:
            * local "up" (direction of zero pitch)
            * local "north" (direction of 90 degree yaw)
            * local "east" (direction of zero yaw)
    :param type: string, either 'pitch' or 'yaw'
    :return: requested angle in degrees
    """

    vector = unit(vector)
    if type == 'pitch':
        angle = safe_acosd(np.vdot(vector, frame[0]))
    elif type == 'yaw':
        inplane = vector - frame[0]*np.vdot(vector, frame[0])
        inplane = unit(inplane)
        angle = safe_acosd(np.vdot(inplane, frame[2]))
        # correct for direction of the angle
        tangential = np.cross(frame[0], frame[2])
        if np.vdot(inplane, tangential) < 0:
            angle = -angle
    else:
        print('Unknown parameter (get_angle_from_frame).\n')
        angle = 0

    if abs(np.imag(angle)) > 0:
        print('-')
    return angle


def safe_acosd(angle):
    """
    Get the arc cosine of the specified angle, in degrees

    :param angle: cosine of the angle (not an angle)
    :return: the angle in degrees
    """
    angle = min(angle,  1)
    angle = max(angle, -1)
    a = np.rad2deg(np.arccos(angle))
    return a

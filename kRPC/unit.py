import numpy as np


def unit(vector):
    """
    Returns a normalized vector. If zero vector is passed, returns this vector.

    :param vector: Any valid XYZ vector.
    :return: Vector of the same direction as v but of magnitude 1.
    """
    if np.linalg.norm(vector) == 0:
        return vector
    else:
        return vector/np.linalg.norm(vector)

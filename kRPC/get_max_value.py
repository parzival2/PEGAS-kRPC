def get_max_value(array):
    """
    Finds (the first occurence of) maximum value in a 1D array.

    :param array: 1D array (there is no safety check on that)
    :return: (index, value)
        index      Index of the maximum element.
        value      Value of the maximum element.
    """
    max_val = max(array)
    idx = 0
    for i in range(len(array)):
        if array[i] == max_val:
            idx = i
            break
    max_i = idx
    return max_i, max_val

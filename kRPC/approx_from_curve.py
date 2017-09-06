def approx_from_curve(x, curve):
    """
    Interpolates a point on a curve given by list of control points.
    If the input is between two points of the curve, output will be
    calculated from the linear function fitted to those points. If the input
    is out of range, extreme values will be output.

    :param x: Argument
    :param curve: Array of 2D (XY) points of shape (n,2)
    :return: Linear approximation for the given argument.
    """

    n = len(curve)

    # If the input is outside the given range, output the extreme value.
    if x > curve[n-1][0]:
        y = curve[n-1][1]
        return y
    elif x < curve[0][0]:
        y = curve[0][1]
        return y

    # Find between which two points defining the curve the argument is.
    for i in range(1, n):
        if curve[i][0] > x:
            break

    # Calculate a linear function coefficients between those points from the
    # given rearrangement:
    #    m*x1+b = y1
    #    m*x2+b = y2
    #    y1-m*x1  =  y2-m*x2
    #    m(x2-x1) = y2-y1
    #    m = (y2-y1)/(x2-x1)
    #    b = y2 - m*x2
    m = (curve[i][1]-curve[i-1][1]) / (curve[i][0]-curve[i-1][0])
    b = curve[i][1] - m*curve[i][0]

    y = m*x+b
    return y

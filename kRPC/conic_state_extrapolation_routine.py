import Global
import numpy as np

mu = Global.mu

# Most of this code is based on the following document
# https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19740018034.pdf


class CSERState:
    # Reference on cser struct:
    #   dtcp        delta tc prim; previous converged value of transfer
    #               time interval
    #   xcp         xc prim, previous converged value of x corresponding to tc
    #   A, D, E     Un function values
    def __init__(self, dtcp, xcp, A, D, E):
        self.dtcp = dtcp
        self.xcp = xcp
        self.A = A
        self.D = D
        self.E = E


def conic_state_extrapolation_routine(r0, v0, dt, last):
    """
    Implementation of Conic State Extrapolation Routine as described by
    Shepperd and Robertson in Space Shuttle GN&C Equation Document 25.
    Starting from the given vehicle state, extrapolates its state into the
    future assuming only central gravity field with no other forces present.
    This implementation has been simulation-proven to give correct results!

    :param r0: Initial position relative to center of the Earth, XYZ
        vector (m).
    :param v0: Initial velocity, XYZ vector (m/s).
    :param dt: Time of extrapolation - finds state so far into the future
        (seconds).
    :param last: Results of previous calculation (needs to be a correct
        struct, can be filled with zeros if this is the first call).
    :return: (r, v, last)
        r    Future position (units like r0).
        v    Future velocity (units like v0).
        last Updated cser struct.
    """

    if last.dtcp == 0:
        dtcp = dt
    else:
        dtcp = last.dtcp

    xcp = last.xcp
    x = xcp        # independent variable for Kepler iteration scheme
    A = last.A
    D = last.D
    E = last.E
    # Program constants
    kmax = 10      # U1 series iterator maximum value
    imax = 10      # Kepler iterator maximum values

    # 5.1 ROUTINE
    # PLATE 5-2 (p24)
    if dt >= 0:
        f0 = 1
    else:
        f0 = -1

    n = 0
    r0m = np.linalg.norm(r0)

    f1 = f0*np.sqrt(r0m/mu)
    f2 = 1/f1
    f3 = f2/r0m
    f4 = f1*r0m
    f5 = f0/np.sqrt(r0m)
    f6 = f0*np.sqrt(r0m)

    ir0 = r0/r0m
    v0s = f1*v0                  # v0 vector with the silly dash on top
    sigma0s = np.vdot(ir0, v0s)  # sigma0 with the silly dash
    b0 = np.vdot(v0s, v0s)-1
    alphas = 1-b0                # alpha with the silly dash

    # PLATE 5-3 (p25)
    xguess = f5*x
    xlast = f5*xcp
    xmin = 0
    dts = f3*dt                  # delta t with the silly dash
    dtlast = f3*dtcp             # delta tlast with the silly dash
    dtmin = 0                    # delta tmin with the silly dash

    # assuming sqrt(alphas) is never smaller than epsilon alpha
    # means orbit is not parabolic (why would it be?)

    xmax = 2*np.pi/np.sqrt(abs(alphas))

    if alphas > 0:
        dtmax = xmax/alphas
        xP = xmax      # xP with the silly dash
        Ps = dtmax     # P with the silly dash
        while dts >= Ps:
            n = n+1
            dts = dts-Ps
            dtlast = dtlast-Ps
            xguess = xguess-xP
            xlast = xlast-xP
    else:
        dtmax, _, _, _ = kepler_transfer_time_interval(xmax, sigma0s, alphas, kmax)
        if dtmax < dts:
            while dtmax >= dts:
                dtmin = dtmax
                xmin = xmax
                xmax = 2*xmax
                dtmax, _, _, _ = kepler_transfer_time_interval(xmax, sigma0s, alphas, kmax)

    # PLATE 5-4 (p26)
    if xmin >= xguess or xguess >= xmax:
        xguess = 0.5*(xmin+xmax)

    dtguess, _, _, _ = kepler_transfer_time_interval(xguess, sigma0s, alphas, kmax)

    if dts < dtguess:
        if xguess < xlast < xmax and dtguess < dtlast < dtmax:
            xmax = xlast
            dtmax = dtlast
    else:
        if xmin < xlast < xguess and dtmin < dtlast < dtguess:
            xmin = xlast
            dtmin = dtlast

    xguess, dtguess, A, D, E = kepler_iteration_loop(imax, dts, xguess, dtguess, xmin, dtmin, xmax, dtmax, sigma0s, alphas, kmax, A, D, E)

    # PLATE 5-5 (p27)
    rs = 1 + 2*(b0*A + sigma0s*D*E)    # r with the silly dash
    b4 = 1/rs

    # The following uses variables (xP, Ps) which are not even created in
    # some execution path (if aplhas>0). However, in the same path variable
    # n is never incremented above 0, so we can restate:
    if n > 0:
        xc = f6*(xguess+n*xP)
        dtc = f4*(dtguess+n*Ps)
    else:
        xc = f6*xguess
        dtc = f4*dtguess

    # Store converged values for the next run
    last.dtcp = dtc
    last.xcp = xc
    last.A = A
    last.D = D
    last.E = E

    # Extrapolated State Vector (ROUTINE 5.3.6, PLATE 5-16 (p38)
    # Implemented inline for simplicity
    F = 1 - 2*A
    Gs = 2*(D*E + sigma0s*A)   # G with the silly dash
    Fts = -2*b4*D*E            # Ft with the silly dash
    Gt = 1 - 2*b4*A

    r = r0m*(F*ir0 + Gs*v0s)
    v = f2*(Fts*ir0 + Gt*v0s)
    return r, v, last


def kepler_transfer_time_interval(xarg, s0s, a, kmax):
    """
    5.3.1 ROUTINE - Kepler Transfer Time Interval
    PLATE 5-9 (p31)

    :param xarg:
    :param s0s:
    :param a:
    :param kmax:
    :return: (t, A, D, E)
        t
        A
        D
        E
    """
    u1 = u1_series_summation(xarg, a, kmax)

    zs = 2*u1                   # z with the silly dash
    E = 1 - 0.5*a*zs**2
    w = np.sqrt(max(0.5+E/2,0)) # added safety check against sqrt from negative value
    D = w*zs
    A = D**2
    B = 2*(E+s0s*D)

    Q = q_continued_fraction(w)

    t = D*(B+A*Q)
    return t, A, D, E


def u1_series_summation(xarg, a, kmax):
    """
    5.3.2 ROUTINE - U1 Series Summation
    PLATE 5-10 (p32)

    :param xarg:
    :param a:
    :param kmax:
    :return: u1
    """
    du1 = 0.25*xarg
    u1 = du1
    f7 = -a*du1**2
    k = 3
    while k < kmax:
        du1 = f7*du1 / (k*(k-1))
        u1old = u1
        u1 = u1+du1
        if u1 == u1old:
            break
        k = k+2
    return u1


def q_continued_fraction(w):
    """
    5.3.3 ROUTINE - Q Continued Fraction
    PLATE 5-11 (p33)

    :param w:
    :return: Q
    """
    if w < 1:
        xq = 21.04-13.04*w
    elif w < 4.625:
        xq = (5/3)*(2*w+5)
    elif w < 13.846:
        xq = (10/7)*(w+12)
    elif w < 44:
        xq = 0.5*(w+60)
    elif w < 100:
        xq = 0.25*(w+164)
    else:
        xq = 70

    # PLATE 5-12 (p34)
    b = 0
    y = (w-1)/(w+1)
    j = np.floor(xq)
    b = y/(1+(j-1)/(j+2)*(1-b))
    while j > 2:
        j = j-1
        b = y/(1+(j-1)/(j+2)*(1-b))

    Q = 1/w**2 * (1 + (2-b/2) / (3*w*(w+1)))
    return Q


def kepler_iteration_loop(imax, dts, xguess, dtguess, xmin, dtmin, xmax, dtmax, s0s, a, kmax, A, D, E):
    """
    5.3.4 ROUTINE - Kepler Iteration Loop
    input skips convergence criteria and "max representable scalar"
    input adds previous values of A, D & E

    :param imax:
    :param dts:
    :param xguess:
    :param dtguess:
    :param xmin:
    :param dtmin:
    :param xmax:
    :param dtmax:
    :param s0s:
    :param a:
    :param kmax:
    :param A:
    :param D:
    :param E:
    :return: (xguess, dtguess, A, D, E)
        xguess
        dtguess
        A
        D
        E
    """
    # PLATE 5-13 (p35)
    i = 1
    while i < imax:
        dterror = dts-dtguess

        if abs(dterror) < 1e-6:  # some arbitrary number, TODO try and look for hints in the paper
            break

        dxs, xmin, dtmin, xmax, dtmax = secant_iterator(dterror, xguess, dtguess, xmin, dtmin, xmax, dtmax)
        xold = xguess
        xguess = xguess + dxs

        if xguess == xold:
            break

        # PLATE 5-14 (p36)
        dtold = dtguess

        dtguess, A, D, E = kepler_transfer_time_interval(xguess, s0s, a, kmax)

        if dtguess == dtold:
            break

        i = i + 1  # this line actually happens on the previous plate
    return xguess, dtguess, A, D, E


def secant_iterator(dterror, xguess, dtguess, xmin, dtmin, xmax, dtmax):
    """
    5.3.5 ROUTINE - Secant Iterator
    skips passing one convergence criterion, instead define it here

    :param dterror:
    :param xguess:
    :param dtguess:
    :param xmin:
    :param dtmin:
    :param xmax:
    :param dtmax:
    :return: (dxs, xmin, dtmin, xmax, dtmax)
        dxs
        xmin
        dtmin
        xmax
        dtmax
    """
    etp = 1e-6
    # PLATE 5-15 (p37)
    dtminp = dtguess-dtmin  # delta tmin prim
    dtmaxp = dtguess-dtmax
    if abs(dtminp) < etp or abs(dtmaxp) < etp:
        dxs = 0
    else:
        if dterror < 0:
            dxs = (xguess-xmax)*(dterror/dtmaxp)
            if (xguess+dxs) <= xmin:
                dxs = (xguess-xmin)*(dterror/dtminp)
            xmax = xguess
            dtmax = dtguess
        else:
            dxs = (xguess-xmin)*(dterror/dtminp)
            if (xguess+dxs) >= xmax:
                dxs = (xguess-xmax)*(dterror/dtmaxp)
            xmin = xguess
            dtmin = dtguess
    return dxs, xmin, dtmin, xmax, dtmax

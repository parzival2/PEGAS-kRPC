import Global

from approx_from_curve import approx_from_curve

g0 = Global.g0


def get_thrust(engines, pressure, time):
    """
    Calculates thrust and closely related parameters for the given set of
    engines at given pressure. Supports engines with constant thrust (eg.
    liquid) and those with thrust profile (eg. solid rocket motors).

    :param engines: Array of struct of type engine.
    :param pressure: Atmospheric pressure (atm).
    :param time: Time since ignition of the engines (seconds).
    :return: (F, dm, isp)
        F          Combined thrust (Newtons).
        dm         Combined mass flow rate (kg/s).
        Isp        Combined specific impulse (seconds).
    """

    n = len(engines)
    p = pressure
    t = time
    F = 0
    dm = 0
    for i in range(n):
        isp1 = engines[i].isp1
        isp0 = engines[i].isp0
        isp = (isp1-isp0)*p+isp0
        if engines[i].mode == 1:
            dm_ = engines[i].flow
        elif engines[i].mode == 2:
            dm_ = approx_from_curve(t, engines[i].data) * engines[i].flow

        dm = dm + dm_
        F = F + isp*dm_*g0
    isp = F/(dm*g0)

    return F, dm, isp

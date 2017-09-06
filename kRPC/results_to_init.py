class InitStruct:
    def __init__(self, t, r, v):
        self.t = t
        self.r = r
        self.v = v
        self.upfg = None


def results_to_init(results):
    """
    Converts a results struct into an initialization struct, allowing easily
    continuing a flight.
    
    :param results: Results struct as returned by flight_sim_3d

    :return: Init struct of type 1.
    """
    n = len(results.plots.t)
    init = InitStruct(results.plots.t[n-1],
                      results.plots.r[n-1],
                      results.plots.v[n-1])
    # Handle UPFG state passing. Copy final state from results, create a
    # dummy state otherwise.
    if results.upfg is not None:
        init.upfg = results.upfg
    else:
        init.upfg = None

    return init

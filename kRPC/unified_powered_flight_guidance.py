import Global
import numpy as np

from conic_state_extrapolation_routine import conic_state_extrapolation_routine, CSERState
from get_angle_from_frame import get_angle_from_frame
from get_thrust import get_thrust
from rodrigues import rodrigues
from unit import unit

g0 = Global.g0

# Most of this code is based on the following document
# https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19740004402.pdf


class Target:
    def __init__(self, angle, normal, radius, velocity):
        self.angle = angle        # angle in radians
        self.normal = normal      # normal of the orbit (opposite angular moment vector)
        self.radius = radius      # magnitude of the radius
        self.velocity = velocity  # magnitude of the velocity


class State:
    def __init__(self, time, mass, radius, velocity):
        self.time = float(time)
        self.mass = float(mass)
        self.radius = radius
        self.velocity = velocity


class UPEGState:
    @classmethod
    def from_state(cls, mu, state, target, theta):
        rd = target.radius * unit(rodrigues(state.radius, -target.normal, theta)) # The idea here is that we want to be in orbit after rotating 20 degrees about the planet
        rgrav = -0.5 * mu * unit(state.radius) / np.linalg.norm(state.radius)**3
        vgo = target.velocity * unit(np.cross(-target.normal, rd)) - state.velocity
        return UPEGState(CSERState(0, 0, 0, 0, 0), np.array([0, 0, 0]), rd, rgrav, 0, state.time, 0, state.velocity, vgo)

    def __init__(self, cser, rbias, rd, rgrav, tb, time, tgo, v, vgo):
        self.cser = cser
        self.rbias = rbias
        self.rd = rd
        self.rgrav = rgrav
        self.tb = tb
        self.time = time
        self.tgo = tgo
        self.v = v
        self.vgo = vgo


class Guidance:
    def __init__(self, pitch, yaw, pitchdot, yawdot, tgo):
        self.pitch = float(pitch)
        self.yaw = float(yaw)
        self.pitchdot = float(pitchdot)
        self.yawdot = float(yawdot)
        self.tgo = float(tgo)


class VehicleStage:
    def __init__(self, mode, m0, engines, maxT, gLim, area, drag):
        self.mode = mode         # thrust mode (1=const thrust, 2=const acc)
        self.m0 = float(m0)      # Initial mass when this stage triggers (includes later stages)
        self.engines = engines   # list of VehicleEngine structures
        self.maxT = float(maxT)  # Max burn time for this stage
        self.gLim = float(gLim)  # g force limit (or 0 for no limit)
        self.area = float(area)  # cross sectional area
        self.drag = drag         # list of drag parameters

    def clone(self):
        engines = [engine.clone() for engine in self.engines]  # deep copy of engines
        drag = None
        if self.drag is not None:  # drag is a list of lists
            drag = [list(item) for item in self.drag]  # deep copy of drag
        return VehicleStage(self.mode, self.m0, engines, self.maxT, self.gLim, self.area, drag)


class VehicleEngine:
    def __init__(self, mode, isp0, isp1, flow, data):
        self.mode = mode         # thrust mode (1=min/max, 2=thrust_profile)
        self.isp0 = float(isp0)  # vacuum isp in seconds
        self.isp1 = float(isp1)  # sea level isp in seconds
        self.flow = float(flow)  # maximum flow rate in kg/s
        self.data = data         # a list containing the min and max throttle settings (mode 1)
                                 # or a list of lists with time and thrust portion (mode 2)

    def clone(self):
        data = None
        if self.data is not None:
            if self.mode == 1:
                data = list(self.data)  # clone the min/max throttle
            elif self.mode == 2:
                data = [list(item) for item in self.data]  # deep copy of thrust profile
        return VehicleEngine(self.mode, self.isp0, self.isp1, self.flow, data)


class DebugState:
    def __init__(self, this, time, r, v, m, dvsensed, vgo1, L1, tgo, L, J, S, Q, P, H, \
                 lambda_vec, rgrav1, rgo1, iz1, rgoxy, rgoz, rgo2, lambdade, lambdadot, \
                 iF, phi, phidot, vthrust, rthrust, vbias, rbias, pitch, EAST, yaw, \
                 rc1, vc1, rc2, vc2, cser_dtcp, cser_xcp, cser_A, cser_D, cser_E, \
                 vgrav, rgrav2, rp, rd, ix, iz2, vd, vgop, dvgo, vgo2, diverge):
        self.this = this
        self.time = time
        self.r = r
        self.v = v
        self.m = m
        self.dvsensed = dvsensed
        self.vgo1 = vgo1
        self.L1 = L1
        self.tgo = tgo
        self.L = L
        self.J = J
        self.S = S
        self.Q = Q
        self.P = P
        self.H = H
        self.lambda_vec = lambda_vec
        self.rgrav1 = rgrav1
        self.rgo1 = rgo1
        self.iz1 = iz1
        self.rgoxy = rgoxy
        self.rgoz = rgoz
        self.rgo2 = rgo2
        self.lambdade = lambdade
        self.lambdadot = lambdadot
        self.iF = iF
        self.phi = phi
        self.phidot = phidot
        self.vthrust = vthrust
        self.rthrust = rthrust
        self.vbias = vbias
        self.rbias = rbias
        self.pitch = pitch
        self.EAST = EAST
        self.yaw = yaw
        self.rc1 = rc1
        self.vc1 = vc1
        self.rc2 = rc2
        self.vc2 = vc2
        self.cser_dtcp = cser_dtcp
        self.cser_xcp = cser_xcp
        self.cser_A = cser_A
        self.cser_D = cser_D
        self.cser_E = cser_E
        self.vgrav = vgrav
        self.rgrav2 = rgrav2
        self.rp = rp
        self.rd = rd
        self.ix = ix
        self.iz2 = iz2
        self.vd = vd
        self.vgop = vgop
        self.dvgo = dvgo
        self.vgo2 = vgo2
        self.diverge = diverge


def unified_powered_flight_guidance(vehicle, target, state, previous):
    """
    Implementation of Unified Powered Flight Guidance in Standard Ascent Mode
    as described by Brand, Brown and Higgins in Space Shuttle GN&C Equation
    Document 24.

    :param vehicle: (Array of) struct defining vehicle performance stage by
        stage. First element of the array should be the currently
        flown stage.
    :param target: Struct defining desired insertion conditions.
    :param state: Struct defining current vehicle physical state.
    :param previous: UPFG internal struct containing results of previous iteration.
    :return: (current, guidance, debug)
        current    UPFG internal struct containing results of this iteration.
        guidance   Struct containing calculated guidance parameters.
        debug      Currently just returns None
    """
    global g0
    
    # Block 0
    gamma	= target.angle     # Desired inertial flight path angle at terminal (cutoff) position
    iy      = target.normal    # Unit vectors relative to desired trajectory plane: iy is normal to desired trajectory plane.
    rdval   = target.radius    # Desired radius magnitude at terminal (cutoff) position
    vdval   = target.velocity  # Desired velocity magnitude at terminal (cutoff) position

    t       = state.time       # Time associated with r, v
    m       = state.mass       # Current estimated vehicle mass
    r       = state.radius     # Vehicle position vector
    v       = state.velocity   # Vehicle velocity vector

    cser    = previous.cser
    rbias   = previous.rbias   # A position bias to account for effects of a rotating thrust vector
    rd      = previous.rd      # Desired terminal (cutoff) position
    rgrav   = previous.rgrav   # Second integral of central force field gravitational acceleration over thrusting maneuver
    tp      = previous.time    # t of previous guidance cycle
    vprev   = previous.v
    vgo     = previous.vgo

    # Block 1
    n  = len(vehicle)   # total number of stages
    SM = [1] * n        # thrust mode (1=const thrust, 2=const acc)
    aL = [0] * n        # acceleration limit for const acceleration mode
    md = [0] * n        # mass flow rate
    ve = [0] * n        # Effective exhaust velocity for phase i
    fT = [0] * n        # thrust
    aT = [0] * n        # acceleration at the beginning of stage
    tu = [0] * n        # "time to burn as if the whole stage was fuel"
    tb = [0] * n        # Estimated burn time remaining in phase i

    for i in range(n):
        SM[i] = vehicle[i].mode
        aL[i] = vehicle[i].gLim * g0
        fT[i], md[i], ve[i] = get_thrust(vehicle[i].engines, 0, 0)
        ve[i] = ve[i] * g0
        aT[i] = fT[i] / vehicle[i].m0
        tu[i] = ve[i] / aT[i]
        tb[i] = vehicle[i].maxT

    # Block 2
    # We need to store dt in order to keep track on maneuver time.
    dt = t - tp  # Guidance cycle time step

    # In the paper, this block assumes the only known thing about vehicle's
    # state change since the last iteration is vector of velocity change
    # (delta-v_sensed) and time interval (delta t). In this implementation
    # however we assume the state is perfectly known and hence the call to
    # Powered Flight Navigation Routine is not necessary.
    # However, we still decrement vgo here!
    dvsensed = v - vprev  # Total velocity change accumulated on accelerometers since last reading
    vgo = vgo - dvsensed  # Velocity to be gained including bias. (vthrust reflects true velocity-to-be-gained)
    vgo1 = vgo

    # Calculation of 'remaining time to burn in stage k' (which here means
    # current stage) is done differently than in the paper. There, t_b,k is
    # decremented by dt every iteration; here this variable is not
    # persistent, but re-read from vehicle data at each iteration, and
    # instead we remember the current burn time tb, so we simply subtract
    # this time from original time-to-burn, obtaining the remaining time.
    tb[0] = tb[0] - previous.tb

    # Block 3
    # Current vehicle parameters have already been obtained in block 1, the
    # only thing different is a_T,k which should be calculated from current
    # mass instead of initial, and subsequently tu,k.
    # This is done according to theory on pages 33-34 and block diagrams on
    # 56-57, although with a small change: original document for the Space
    # Shuttle assumed that standard ascent will be finalized with a
    # predetermined OMS burn (Orbiter's SSMEs burning fuel from ET will only
    # take it so far, then the ET is jettisoned and vehicle coasts for a
    # predetermined amount of time (tc), after which the last burn phase
    # circularizes). Therefore, active guidance did not calculate the
    # integral Li(n). Stages 1..n-2 were assumed to burn out completely
    # (hence the logarithmic expression for them), and stage n-1 has burn
    # time set accordingly, so the total delta-v expended equals vgo.
    # In this application however, there is no such thing as a predetermined
    # final stage velocity. Therefore, stages n-1 are assumed to burn out
    # completely, and the last one is adjusted, so it burns out only as long
    # as needed.

    if SM[0] == 1:
        aT[0] = fT[0] / m
    elif SM[0] == 2:
        aT[0] = aL[0]

    tu[0] = ve[0] / aT[0]
    L = 0
    Li = [0]*n
    for i in range(n-1):
        if SM[i] == 1:
            Li[i] = ve[i]*np.log(tu[i] / (tu[i]-tb[i]))
        elif SM[i] == 2:
            Li[i] = aL[i]*tb[i]

        L = L + Li[i]
        # If we have more stages than we need to get to orbit, redo the
        # whole calculation but skip the last stage.
        if L > np.linalg.norm(vgo):
            return unified_powered_flight_guidance(vehicle[0:n-1], target, state, previous)
    Li[n-1] = np.linalg.norm(vgo) - L
    # Now for each stage its remaining time of burn is calculated (tbi) and
    # in the same pass accumulated into a total time-to-go of the maneuver.
    tgoi = [0]*n  # Time-to-go until end of ith phase
    for i in range(n):
        if SM[i] == 1:
            tb[i] = tu[i]*(1-np.exp(-Li[i]/ve[i]))
        elif SM[i] == 2:
            tb[i] = Li[i] / aL[i]
        if i == 0:
            tgoi[i] = tb[i]
        else:
            tgoi[i] = tgoi[i-1] + tb[i]
    L1 = Li[0]
    tgo = tgoi[n-1]

    # Block 4
    L = 0
    J = 0
    Ji = [0]*n
    S = 0
    Si = [0]*n
    Q = 0
    Qi = [0]*n
    H = 0
    P = 0
    Pi = [0]*n
    # Major loop of the whole block, almost exactly as in the block diagrams.
    for i in range(n):
        # Variable tgoi1 represents t_go,i-1 only is determined in a safe
        # way (as to not exceed the array).
        if i == 0:
            tgoi1 = 0
        else:
            tgoi1 = tgoi[i-1]

        # Constant thrust vs constant acceleration mode
        if SM[i] == 1:
            Ji[i] = tu[i]*Li[i] - ve[i]*tb[i]
            Si[i] = -Ji[i] + tb[i]*Li[i]
            Qi[i] = Si[i]*(tu[i]+tgoi1) - 0.5*ve[i]*tb[i]**2
            Pi[i] = Qi[i]*(tu[i]+tgoi1) - 0.5*ve[i]*tb[i]**2 * ((1.0/3)*tb[i]+tgoi1)
        elif SM[i] == 2:
            Ji[i] = 0.5*Li[i]*tb[i]
            Si[i] = Ji[i]
            Qi[i] = Si[i]*((1.0/3)*tb[i] + tgoi1)
            Pi[i] = (1.0/6)*Si[i]*(tgoi[i]**2 + 2*tgoi[i]*tgoi1 + 3*tgoi1**2)

        # Common for both modes
        Ji[i] = Ji[i] + Li[i]*tgoi1
        Si[i] = Si[i] + L*tb[i]
        Qi[i] = Qi[i] + J*tb[i]
        Pi[i] = Pi[i] + H*tb[i]

        # No coast period before the last stage.

        L = L + Li[i]
        J = J + Ji[i]
        S = S + Si[i]
        Q = Q + Qi[i]
        P = P + Pi[i]
        H = J*tgoi[i] - Q

    # Block 5
    lambda_vec = unit(vgo)     # Unit vector in direction of vgo
    # print('lambda %s; vgo = %s' % (lambda_vec, vgo))
    rgrav1 = rgrav
    if not np.isclose(previous.tgo, 0):
       rgrav = (tgo/previous.tgo)**2 * rgrav
    rgo = rd - (r + v*tgo + rgrav)
    rgo1 = rgo
    iz = unit(np.cross(rd, iy))
    # print('iz = %s; rd = %s; iy = %s' % (iz, rd, iy))
    iz1 = iz
    rgoxy = rgo - np.vdot(iz, rgo)*iz
    rgoz = (S - np.vdot(lambda_vec, rgoxy)) / np.vdot(lambda_vec, iz)
    rgo = rgoxy + rgoz*iz + rbias
    lambdade = Q - S*J/L
    lambdadot = (rgo - S*lambda_vec) / lambdade  # Time derivative of unit vector coincident with lambda_vec, but rotating with desired thrust vector turning rate omega_f
    iF = unit(lambda_vec - lambdadot*J/L)     # Unit thrust vector
    phi = np.arccos(np.vdot(iF, lambda_vec))
    phidot = -phi*L/J
    vthrust = (L - 0.5*L*phi**2 - J*phi*phidot - 0.5*H*phidot**2)*lambda_vec  # First integral of thrust acceleration vector over thrusting maneuver
    vthrust = vthrust - (L*phi + J*phidot)*unit(lambdadot)
    # print("vthrust = %s" % vthrust)
    rthrust = (S - 0.5*S*phi**2 - Q*phi*phidot - 0.5*P*phidot**2)*lambda_vec  # Second integral of thrust acceleration vector over thrusting maneuver
    rthrust = rthrust - (S*phi + Q*phidot)*unit(lambdadot)
    vbias = vgo - vthrust  # A velocity bias to account for the effects of a rotating thrust vector (vbias = vgo - vthrust)
    rbias = rgo - rthrust  # A position bias to account for the effects of a rotating thrust vector (rbias = rgo - rthrust)

    # Block 6 - original document does not contain any implementation
    # TODO - pitch and yaw RATES
    UP = unit(r)
    NORTH = np.array([0, 0, 1])
    EAST = unit(np.cross(NORTH, UP))
    frame = [UP, NORTH, EAST]
    # print('Frame: %s; iF = %s' % (frame, iF))
    pitch = get_angle_from_frame(iF, frame, 'pitch')
    yaw = get_angle_from_frame(iF, frame, 'yaw')

    # Block 7 - this calls the Conic State Extrapolation Routine
    rc1 = r - 0.1*rthrust - (1.0/30)*vthrust*tgo  # Vehicle position vector at beginning of gravity computation coast segment
    vc1 = v + 1.2*rthrust/tgo - 0.1*vthrust     # Vehicle velocity vector at beginning of gravity computation coast segment
    # Vehicle velocity vector at end of gravity computation coast segment
    # Vehicle position vector at end of gravity computation coast segment
    rc2, vc2, cser = conic_state_extrapolation_routine(rc1, vc1, tgo, cser)
    vgrav = vc2 - vc1            # First integral of central force field gravity acceleration over thrusting maneuver
    rgrav = rc2 - rc1 - vc1*tgo  # Second integral of central force field gravitational acceleration over thrusting maneuver

    # Block 8
    rho = 0  # Damping factor used in determining the change in dvgo - see 4-20 for a proper definition
    rp = r + v*tgo + rgrav + rthrust
    rp = rp - np.vdot(rp, iy)*iy
    rd = rdval*unit(rp)
    ix = unit(rd)
    iz = np.cross(ix, iy)
    vd = vdval*np.matmul(np.transpose([ix, iy, iz]), [np.sin(gamma), 0, np.cos(gamma)])
    vgop = vd - v - vgrav + vbias
    dvgo = rho*(vgop-vgo)  # big values (0.8+) cause bananas; standard ascent uses 0 (?)
    # print('vd = %s; gamma = %f; vgop = %s; v = %s; vgrav = %s; vbias = %s' % (vd, gamma, vgop, v, vgrav, vbias))
    vgo = vgop + dvgo  # 5-15 shows vgo + dvgo, but 4-21 shows vgop + dvgo

    current = previous
    current.cser = cser
    current.rbias = rbias
    current.rd = rd
    current.rgrav = rgrav
    current.tb = current.tb + dt
    current.time = t
    current.tgo = tgo
    current.v = v
    current.vgo = vgo

    guidance = Guidance(pitch, yaw, 0, 0, tgo)

    debug = DebugState(0, t, r, v, m, dvsensed, vgo1, L1, tgo, L, J, S, Q, P, H,
                       lambda_vec, rgrav1, rgo1, iz1, rgoxy, rgoz, rgo, lambdade, lambdadot,
                       iF, phi, phidot, vthrust, rthrust, vbias, rbias, pitch, EAST, yaw,
                       rc1, vc1, rc2, vc2, cser.dtcp, cser.xcp, cser.A, cser.D, cser.E,
                       vgrav, rgrav, rp, rd, ix, iz, vd, vgop, dvgo, vgo, 0)

    return current, guidance, debug

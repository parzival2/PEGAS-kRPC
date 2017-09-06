import numpy as np
import time
import krpc
import init_simulation
import unified_powered_flight_guidance as upfg
from launch_targeting import launch_targeting, LaunchSite
from flight_manager import flight_manager
from flight_sim_3d import GravityTurnControl, cart2sph

g0 = init_simulation.g0
mu = init_simulation.mu
R = init_simulation.R
period = init_simulation.period

conn = None
space_center = None
vessel = None


def init_from_ksp():
    global conn
    global space_center
    global vessel
    global mu
    global R
    global period
    conn = krpc.connect(name='Launch to orbit')
    space_center = conn.space_center
    vessel = space_center.active_vessel
    mu = vessel.orbit.body.gravitational_parameter
    R = vessel.orbit.body.equatorial_radius
    period = vessel.orbit.body.rotational_period


def vec_yz(v):
    return np.array([v[0], v[2], v[1]])


def get_engines(vessel, decouple_stage):
    global g0
    engines = []
    for engine in vessel.parts.engines:
        if engine.part.decouple_stage == decouple_stage:
            flow = engine.max_vacuum_thrust / (engine.vacuum_specific_impulse * g0)
            # TODO: determine throttle limits on the engine
            engines.append(upfg.VehicleEngine(1, engine.vacuum_specific_impulse, engine.kerbin_sea_level_specific_impulse, flow, [0, 1]))
    return engines


def get_max_stage(part_list):
    max_stage = 0
    for part in part_list:
        if part.stage > max_stage:
            max_stage = part.stage
    return max_stage


def get_mass(vessel, decouple_stage):
    mass = 0
    for part in vessel.parts.all:
        if part.decouple_stage == decouple_stage:
            mass = mass + part.mass
    return mass


def get_fuel_mass(vessel, decouple_stage):
    fuel_mass = 0
    for part in vessel.parts.all:
        if part.decouple_stage == decouple_stage:
            fuel_mass = fuel_mass + part.mass - part.dry_mass
    return fuel_mass


def check_engine(vessel):
    """
    Check to see if there are active engines

    :param vessel: Vessel object from kRPC
    :return: whether any engines are active and have thrust available
    """
    for engine in vessel.parts.engines:
        if engine.active and engine.available_thrust > 0:
            return True
    return False


def stage_if_needed(vessel):
    """
    Check to see if we need to stage, and if so, activate the next stage

    :param vessel: Vessel object from kRPC
    :return: whether we activated a new stage
    """
    if check_engine(vessel):
        return False

    print('There is no active engine, checking Propellant condition')
    for engine in vessel.parts.engines:
        if engine.has_fuel:
            for engine_module in engine.part.modules:
                if engine_module.has_field('Propellant'):
                    if engine_module.get_field('Propellant') == 'Very Stable':
                        print('Engine is ready; staging')
                        vessel.control.forward = 0
                        vessel.control.throttle = 1
                        vessel.control.activate_next_stage()
                        return True

    print('No engine is ready')
    return False


def apply_guidance(pitch, yaw, desired_throttle):
    global vessel
    # KSP target_pitch has 90 as up and 0 as flat, while guidance.pitch has 0 as up and 90 as flat
    vessel.auto_pilot.target_pitch = 90 - pitch
    vessel.auto_pilot.target_heading = 90 - yaw
    vessel.control.throttle = desired_throttle
    time.sleep(0.2)  # run our guidance 5x a second
    # See if we need to auto-stage
    stage_if_needed(vessel)
    # we do this after the sleep, to get the last tick of thrust
    vessel.control.throttle = desired_throttle
    # print('pitch: %f, yaw: %f, throttle: %f' % (pitch, yaw, desired_throttle))


def get_state():
    global liftoff_time
    global vessel
    r = vec_yz(vessel.flight(vessel.orbit.body.non_rotating_reference_frame).center_of_mass)
    v = vec_yz(vessel.flight(vessel.orbit.body.non_rotating_reference_frame).velocity)
    m = vessel.mass
    t = vessel.met - liftoff_time
    # print('current state: r=%s, v=%s, m=%f, t=%f' % (r, v, m, t))
    return r, v, m, t


def site_from_position(vessel):
    # The initial site is based on the rotating reference frame, but we need a position in the non-rotating reference frame
    r = vec_yz(vessel.flight(vessel.orbit.body.non_rotating_reference_frame).center_of_mass)
    longitude, latitude, r = cart2sph(r[0], r[1], r[2])
    updated_site = LaunchSite(np.rad2deg(latitude), np.rad2deg(longitude), r - vessel.orbit.body.equatorial_radius)
    return updated_site


def warp_to_launch_time(space_center, launch_time):
    print('Launch time is %f' % launch_time)
    game_launch_time = space_center.ut + launch_time
    space_center.warp_to(game_launch_time - 10)

    while (space_center.ut - game_launch_time) < 0:
        print('Time to launch %f' % (space_center.ut - game_launch_time))
        time.sleep(1)


def vehicle_from_vessel(vessel, gravity_turn_t):
    """
    Generate a basic vehicle structure from the provided vessel.
    More complicated flights will require manual construction of this structure.

    :param vessel: kRPC Vessel object
    :param gravity_turn_t: Time to spend in the initial gravity turn (s)
    :return: list of VehicleStage objects representing the various flight stages
    """
    vehicle_stages = []
    max_stage = get_max_stage(vessel.parts.all)
    engines = [None]*(max_stage+2)
    masses = [None]*(max_stage+2)
    fuel_masses = [None]*(max_stage+2)
    for decouple_stage in range(-1, max_stage+1):
        engines[decouple_stage+1] = get_engines(vessel, decouple_stage)
        masses[decouple_stage+1] = get_mass(vessel, decouple_stage)
        fuel_masses[decouple_stage+1] = get_fuel_mass(vessel, decouple_stage)
    # Assumption here is that we jettison the fuel tank at the same times as their corresponding engines
    for stage in range(max_stage+1, -1, -1):
        if len(engines[stage]) == 0:
            continue
        fuel_mass = fuel_masses[stage]
        flow = reduce(lambda x, engine: x + engine.flow, engines[stage], 0.0)
        stage_time = fuel_mass / flow
        mass = 0.0
        for i in range(0, stage+1):
            mass = mass + masses[i]
        vehicle_stages.append(upfg.VehicleStage(1, mass, engines[stage], stage_time, 0, 0, [[0, 0]]))
    if gravity_turn_t > 0.0:
        split_stage = vehicle_stages[0]
        split_stage_flow = reduce(lambda x, engine: x + engine.flow, split_stage.engines, 0.0)
        gravity_stage = upfg.VehicleStage(1, split_stage.m0, split_stage.engines, gravity_turn_t,
                                          split_stage.gLim, split_stage.area, split_stage.drag)
        upfg_stage = upfg.VehicleStage(1, split_stage.m0 - split_stage_flow * gravity_turn_t,
                                       split_stage.engines, split_stage.maxT - gravity_turn_t,
                                       split_stage.gLim, split_stage.area, split_stage.drag)
        vehicle_stages = [gravity_stage, upfg_stage] + vehicle_stages[1:]
    return vehicle_stages


# Connect to ksp, and set up some globals
init_from_ksp()

# Generate a basic vehicle description, with a 90 second gravity turn, followed by UPFG for each stage
vehicle = vehicle_from_vessel(vessel, 90)

# Target orbit parameters
apoapsis = 250
periapsis = 170
inclination = 28.61  # 28.61 = 23.45+5.16
lan = 0  # Set this to None if you don't care

print('Launching into orbit with apoapsis of %fkm, periapsis %fkm, inclination %f and LAN %s degrees' %
      (apoapsis, periapsis, inclination, lan))

launch_site = LaunchSite(vessel.flight().latitude, vessel.flight().longitude, 0)
if inclination < launch_site.latitude:
    print('Target inclination is below site latitude. Setting inclination to match latitude.')
    inclination = launch_site.latitude

# Determine now long we need to wait before launching, to make our longitude of ascending node match
if lan is not None:
    estimated_lan, _, _ = launch_targeting(launch_site, periapsis, apoapsis, inclination, 2)
    # adjust lan based on the current rotation of the celestial body
    estimated_lan = np.rad2deg(vessel.orbit.body.rotation_angle) + estimated_lan
    wait_time = period * np.mod(lan - estimated_lan + 360, 360) / 360
    warp_to_launch_time(space_center, wait_time)

liftoff_time = vessel.met
vessel.control.throttle = 1

# If we don't have any engines going, stage until we do
while not check_engine(vessel):
    stage_if_needed(vessel)
    time.sleep(0.1)

# Wait until we have enough thrust to lift off
while vessel.thrust < vessel.mass * mu / vessel.orbit.radius**2:
    time.sleep(0.01)

# Release launch clamps
while len(vessel.parts.launch_clamps) > 0:
    print('Releasing launch clamp')
    vessel.control.activate_next_stage()
    time.sleep(0.1)

# Wait until our thrust is at full, before we start maneuvers
while vessel.thrust < vessel.available_thrust:
    time.sleep(0.01)

vessel.auto_pilot.engage()
vessel.auto_pilot.attenuation_angle = (1, 1, 0.4)
vessel.auto_pilot.target_heading = 0
vessel.auto_pilot.target_pitch = 90

print('Proceeding Launch..')

while vessel.flight(vessel.orbit.body.reference_frame).speed < 30:
    time.sleep(0.1)

print('Clear from launch tower..')
print('Begin Pitch and Roll Program..')
vessel.auto_pilot.target_roll = 0

# Recompute the target information, now that we've advanced time, and launched
site = site_from_position(vessel)
_, azm, target = launch_targeting(site, periapsis, apoapsis, inclination, 2.0)

stage1 = GravityTurnControl(8.5, 50, azm)

STS = flight_manager(vehicle, site, target, .2, stage1, 2, 0, [], apply_guidance, get_state)
vessel.control.throttle = 0.0

print('Mission Success')

vessel.auto_pilot.disengage()
conn.close()


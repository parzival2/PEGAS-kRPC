import unified_powered_flight_guidance as upfg
import space_shuttle_thrust_profile
import init_simulation
from vehicle_tools import vehicle_tools

g0 = init_simulation.g0

# Full model of the Space Shuttle, 4 stages:
#    SSME + SRB
#    constant-thrust SSME mode
#    acceleration-limiting SSME mode
#    OMS
# Data taken from:
#    http://www.braeunig.us/space/specs/shuttle.htm
#    http://www.braeunig.us/space/specs/orbiter.htm
#    http://www.astronautix.com/s/srb.html
#    https://en.wikipedia.org/wiki/Space_Shuttle_orbiter#Shuttle_Orbiter_Specifications_.28OV-105.29
orbiter = 79135        # OV-105 Endeavour, full OMS/RCS fuel tanks
payload = 25000        # max for Endeavour

thrust_profile = space_shuttle_thrust_profile.thrust_profile

# The Space Shuttle is modeled with four stages
vehicle = [None]*4

# SRB+SSME
stage_m0 = 2*587000 + 765000 + orbiter + payload   # launch mass [kg]   | SRBs + ET (35t structure, 730t fuel) + orbiter + payload
stage_engines = [None]*2
stage_engines[0] = upfg.VehicleEngine(1, 452, 366, 1462.7, [0.670, 1.045])  # 3 SSMEs
stage_engines[1] = upfg.VehicleEngine(2, 269, 237, 2*13500000/(242*g0),   # approx. peak thrust from STS-107 mission plots
                                      thrust_profile)  # 2 SRBs
stage_time = 124                                       # SRB burn time [s]
stage_area = 2*10.8 + 55.4 + 15.2                      # cross section [m2]
stage_drag = [[0.0, 0.08],
              [250, 0.08],
              [343, 1.20],
              [999, 0.50],
              [9999, 0.40]]                        # drag curve - not supported by any real data!
stage = upfg.VehicleStage(1, stage_m0, stage_engines, stage_time, 0, stage_area, stage_drag) # constant thrust
vehicle[0] = stage.clone()

# SSME const-thrust
stage_m0 = 765000 + orbiter + payload
stage.engines = stage.engines[:-1]
stage_engines = stage_engines[:-1]                                # SRBs were jettisoned
stage_m0 = stage_m0 - vehicle_tools('mass', stage, 1, stage.maxT) # SSMEs burned some of the fuel in ET away
stage_time = 320                                                  # acceleration exceeds 3G after that time
stage_area = 55.4 + 15.2
stage = upfg.VehicleStage(1, stage_m0, stage_engines, stage_time, 0, stage_area, stage_drag)
vehicle[1] = stage.clone()

# SSME g-limited
stage_m0 = stage_m0 - vehicle_tools('mass', stage, 1, stage.maxT)
stage_fuel = 730000 - vehicle_tools('mass', stage, 1, 124+320)  # fuel left after previous stages
stage = upfg.VehicleStage(2, stage_m0, stage_engines, 0, 3, stage_area, stage_drag)  # constant acceleration; acceleration limit of 3 Gs
stage.maxT = vehicle_tools('tgo', stage, 2, [stage_fuel, stage.gLim])
vehicle[2] = stage.clone()

# OMS
stage_m0 = orbiter + payload
stage_area = 15.2
stage_engines = [None]
stage_engines[0] = upfg.VehicleEngine(1, 313, 313, 2*26700/(313*g0), thrust_profile)
stage = upfg.VehicleStage(1, stage_m0, stage_engines, 0, 0, stage_area, stage_drag)
stage.maxT = vehicle_tools('tgo', stage, 1, 8174+13486) # mass of MMH and N2O4 in OMS/RCS pods
vehicle[3] = stage.clone()

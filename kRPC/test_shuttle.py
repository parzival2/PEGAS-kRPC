from launch_targeting import launch_targeting
from space_shuttle import vehicle
from flight_manager import flight_manager
from create_launch_site import create_launch_site
from flight_sim_3d import GravityTurnControl


site = create_launch_site('Kennedy')

periapsis = 100
apoapsis = 250
inclination = 51.65
lan, azm, target = launch_targeting(site, periapsis, apoapsis, inclination, 2.0)

stage1 = GravityTurnControl(8.5, 50, azm)

STS = flight_manager(vehicle, site, target, 0.2, stage1, 2, 0, [])

# telemetry(STS.powered, STS.coast, 1)
# dbgIntegrals(STS.powered(2:STS.n), 2)
# trajectory(STS.powered, STS.coast, target, 2, 1, 3)

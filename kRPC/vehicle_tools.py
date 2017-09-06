import numpy as np


def vehicle_tools(task, stage, mode, input):
    """
    VEHICLETOOLS calculate often-used parameters for a given stage and mode
    (constant thrust of constant acceleration). Exact input depends on type of
    call. Constant acceleration mode can only be used for single-engine stages
    (collapse multiple engines with the same parameters into one) that are
    throttleable (ie. no thrust profile allowed for constant acceleration).
    The function will not check for violation of throttle limits (will warn if
    they are exceeded though). Constant acceleration assumes vacuum conditions
    for Isp calculation.
    
    :param stage: Struct of type vehicle stage (ie. an element of a vehicle
        array).
    :param mode: Mode to be used; allowed choices:
        1 - constant thrust
        2 - constant acceleration
    
    time = vehicle_tools('tgo', stage, mode, fuel)
    # Calculates maximum burn time for a stage with given set of engines and
    # total mass of fuel.

%INPUT
%    fuel       If mode==1: Mass of fuel to be burned (kg).
%               If mode==2: Array of length 2, containing the following:
%                               fuel(1) = mass of fuel to be burned (kg).
%                               fuel(2) = acceleration limit (G).
%
%OUTPUT
%    time       Time it takes to burn the given amount of fuel.
%
%
%mass = VEHICLETOOLS('mass', stage, mode, time)
%Calculates mass of fuel burned by a stage with given set of engines over
%some given time.
%
%INPUT
%    time       If mode==1: Time of burn (s).
%               If mode==2: Array of length 2, containing the following:
%                               time(1) = time of burn (s).
%                               time(2) = acceleration limit (G).
%OUTPUT
%    mass       Mass of fuel burned over the given time.
    """
    if task == 'tgo':
        if mode == 1:  # constant-thrust
            dm = 0
            for engine in stage.engines:
                if engine.mode != 1:
                    print('vehicle_tools couldnt calculate tgo for the given stage. One of the engines is in profiled or unknown mode.')
                    return 0
                dm = dm + engine.flow
            return input / dm
        elif mode == 2: # constant-acceleration
            if len(stage.engines) != 1:
                print('vehicle_tools couldnt calculate tgo for the given stage. Only single-engined stages can run in constant-acceleration mode.')
                return 0
            m0 = stage.m0
            isp = stage.engines[0].isp0
            result = isp / input[1] * np.log(m0/(m0-input[0]))
            # check for violation of engine throttling limits
            max_flow = input[1] * m0 / isp
            min_flow = input[1] * (m0-input[0]) / isp
            nominal_flow = stage.engines[0].flow
            min_throttle = stage.engines[0].data[0]
            max_throttle = stage.engines[0].data[1]
            if max_flow > max_throttle * nominal_flow or min_flow > max_throttle * nominal_flow:
                print('vehicle_tools found a violation of upper throttle limit! Time calculated by this call will likely underestimate true burn time.')
            if min_flow < min_throttle * nominal_flow or max_flow < min_throttle * nominal_flow:
                print('vehicle_tools found a violation of lower throttle limit! Time calculated by this call will likely overestimate true burn time.')
            return result
        else:
            print('vehicle_tools couldnt calculate tgo for the given stage. Unknown throttle mode.')
            return 0
    elif task == 'mass':
        if mode == 1:    # constant-thrust
            dm = 0
            for engine in stage.engines:
                if engine.mode != 1:
                    print('vehicle_tools couldn\'t calculate burned mass for the given stage. One of the engines is in profiled or unknown mode.')
                    return 0
                dm = dm + engine.flow
            return input * dm
        elif mode == 2:  #constant-acceleration
            if len(stage.engines) != 1:
                print('vehicle_tools couldn\'t calculate burned mass for the given stage. Only single-engined stages can run in constant-acceleration mode.')
                return 0
            m0 = stage.m0
            isp = stage.engines[0].isp0
            alim = input[1]
            t = input[0]
            return m0 - m0 * np.exp(-t / (isp/alim))
        else:
            print('vehicle_tools couldn\'t calculate burned mass for the given stage. Unknown throttle mode.')
            return 0
    else:
        print('vehicle_tools failed to understand the task.')
        return 0

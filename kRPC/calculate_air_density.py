def calculate_air_density(pressure, temperature):
    """
    Calculates air density from ideal gas law.
    
    :param pressure: Air pressure in Pascals
    :param temp: Air temperature in Kelvins
    :return density: Air density in kg/m^2
    """

    R = 286.9  # specific gas constant for air [J/(kg*K)]
    density = pressure/(R*temperature)
    return density

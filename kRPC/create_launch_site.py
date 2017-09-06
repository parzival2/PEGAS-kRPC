from launch_targeting import LaunchSite


def create_launch_site(name):
    """
    site = CREATELAUNCHSITE(name)
    Outputs a properly formatted init struct from a given string representing
    a launch site name.

   :param name: Short name of a launch site. Supported choices:
                full name                           short name
                ------------------------------------------------
                Kennedy Space Center                KSC, Kennedy
                Vandenberg Air Force Base           Vandenberg
                Guiana Space Centre                 Kourou
                Baikonur Cosmodrome                 Baikonur
                Plesetsk Cosmodrome                 Plesetsk
                Uchinoura Space Center              Uchinoura
                RAAF Woomera Test Range             Woomera
                Jiuquan Satellite Launch Center     Jiuquan
 
    :return: 
     site       Init struct of type 0 containing longitude, latitude and
                altitude of a selected launch site.
    """

    if name == 'KSC':
        ls = LaunchSite(28.52406, -80.65085, 0)
    elif name == 'Kennedy':
        ls = create_launch_site('KSC')
    elif name == 'Vandenberg':
        ls = LaunchSite(34.75083, -120.49778, 0)
    elif name == 'Kourou':
        ls = LaunchSite(5.15972, -52.65028, 0)
    elif name == 'Baikonur':
        ls = LaunchSite(45.91194, 63.31028, 92)
    elif name == 'Plesetsk':
        ls = LaunchSite(62.92556, 40.57778, 131)
    elif name == 'Uchinoura':
        ls = LaunchSite(31.25194, 131.08914, 0)
    elif name == 'Woomera':
        ls = LaunchSite(-30.94907, 136.53418, 137)
    elif name == 'Jiuquan':
        ls = LaunchSite(41.11803, 100.46330, 1045)
    else:
        ls = create_launch_site('KSC')
        print('Unknown launch site - returning KSC.')
    return ls

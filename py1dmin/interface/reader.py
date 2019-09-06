"""
reads info
"""

import autoparse.pattern as app
import autoparse.find as apf


def lennard_jones(output_string):
    """ reads the lennard jones params from the output
    """

    sigma_ptt = ('SIGMA' + app.SPACES + '=' + app.SPACES +
                 app.capturing(app.FLOAT))
    epsilon_ptt = ('EPSILON' + app.SPACES + '=' + app.SPACES +
                   app.capturing(app.FLOAT))

    sigma = float(apf.last_capture(sigma_ptt, output_string))
    epsilon = float(apf.last_capture(epsilon_ptt, output_string))

    return sigma, epsilon


def combine_params(param1, param2, rule='default'):
    """ perform a combining rule for two parameters
    """

    if rule == 'default':
        combined_param = (param1 + param2) / 2.0
    else:
        raise NotImplementedError

    return combined_param

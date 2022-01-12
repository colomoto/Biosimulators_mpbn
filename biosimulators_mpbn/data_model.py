""" Data model for the mapping between KISAO terms and simulation methods and their arguments

:Author: Loïc Paulevé <loic.pauleve@labri.fr>
:Date: 2022-01-12
:Copyright: 2022, Loïc Paulevé
:License: MIT
"""

from biosimulators_utils.sedml.data_model import SteadyStateSimulation
import collections
import biosimulators_mpbn

__all__ = ['KISAO_METHOD_MAP']

KISAO_METHOD_MAP = collections.OrderedDict([
    ('KISAO_0000660', { # TODO create specific term
        'kisao_id': 'KISAO_0000660',
        'name': 'logical model stable state search method',
        'simulation_type': SteadyStateSimulation,
        'method': 'fixedpoints',
        'method_args': lambda simulation: [],
        'parameters': {},
    }),
    ('KISAO_0000661', { # TODO create specific term
        'kisao_id': 'KISAO_0000661',
        'name': 'logical model trap space identification method',
        'simulation_type': SteadyStateSimulation,
        'method': 'attractors',
        'method_args': lambda simulation: [],
        'parameters': {},
    }),
])

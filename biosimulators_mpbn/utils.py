
import os
from biosimulators_utils.report.data_model import VariableResults
from biosimulators_utils.sedml.data_model import (  # noqa: F401
    ModelLanguage, )
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.warnings import warn, BioSimulatorsWarning
import biosimulators_utils.model_lang.sbml.utils
import biosimulators_utils.sedml.validation
import biosimulators_utils.xml.utils
from kisao.utils import get_preferred_substitute_algorithm_by_ids

from .data_model import KISAO_METHOD_MAP  # noqa: F401

import lxml.etree
import numpy

from ginsim.gateway import japi as ginsim_japi
import biolqm
import mpbn

__all__ = [
    "read_model",
    'get_variable_target_xpath_ids',
    'get_variable_results',
    'exec_simulation',
    'set_up_simulation',
    'validate_variables',
]

def read_model(filename, language):
    """ Read a model

    Args:
        language (:obj:`ModelLanguage`): language

    Returns:
        :obj:`py4j.java_gateway.JavaObject`: model
    """
    if language == ModelLanguage.SBML:
        format = 'sbml'
    else:
        format = None

    if not os.path.isfile(filename):
        raise FileNotFoundError('`{}` is not a file.'.format(filename))

    model = ginsim_japi.lqm.load(filename, format)

    if model is None:
        raise ValueError('Model `{}` could not be loaded.'.format(filename))

    model = mpbn.MPBooleanNetwork(biolqm.to_minibn(model))

    return model

def set_up_simulation(simulation, config=None):
    """ Set up an analysis

    Args:
        simulation (:obj:`Simulation`): analysis
        config (:obj:`Config`, optional): configuration

    Returns:
        :obj:`tuple`:

            * :obj:`str`: KiSAO of algorithm to execute
            * :obj:`str`: name of the :obj:`biolqm` simulation/analysis method
            * :obj:`types.LambdaType` of :obj:`Simulation` -> :obj:`list` of :obj:`str`: arguments for simulation method
    """
    # simulation algorithm
    alg_kisao_id = simulation.algorithm.kisao_id
    alg_substitution_policy = get_algorithm_substitution_policy(config=config)
    exec_kisao_id = get_preferred_substitute_algorithm_by_ids(
        alg_kisao_id, KISAO_METHOD_MAP.keys(),
        substitution_policy=alg_substitution_policy)
    alg_props = KISAO_METHOD_MAP[exec_kisao_id]
    if simulation.__class__ != alg_props['simulation_type']:
        raise NotImplementedError('{} simulations cannot be executed with `{}` ({}).'.format(
            simulation.__class__.__name__, exec_kisao_id, alg_props['name']))

    method = alg_props['method']
    method_args = []

    # Apply the algorithm parameter changes specified by `simulation.algorithm.parameter_changes`
    if exec_kisao_id == alg_kisao_id:
        for change in simulation.algorithm.changes:
            param_props = alg_props['parameters'].get(change.kisao_id, None)

            if param_props is None:
                if (
                    ALGORITHM_SUBSTITUTION_POLICY_LEVELS[alg_substitution_policy]
                    > ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
                ):
                    warn('Unsuported algorithm parameter `{}` was ignored.'.format(change.kisao_id), BioSimulatorsWarning)
                else:
                    raise NotImplementedError('Algorithm parameter `{}` is not supported.'.format(change.kisao_id))
            else:
                if validate_str_value(change.new_value, param_props['type']):
                    parsed_value = parse_value(change.new_value, param_props['type'])
                    method_args.extend(param_props['method_args'](parsed_value))

                else:
                    if (
                        ALGORITHM_SUBSTITUTION_POLICY_LEVELS[alg_substitution_policy]
                        > ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
                    ):
                        msg = 'Unsuported algorithm parameter value `{}` of `{}` was ignored. The value must be a {}.'.format(
                            change.new_value, change.kisao_id, param_props['type'].value)
                        warn(msg, BioSimulatorsWarning)
                    else:
                        msg = '`{}` is not a valid value of algorithm parameter `{}`. The value must be a {}.'.format(
                            change.new_value, change.kisao_id, param_props['type'].value)
                        raise ValueError(msg)
    else:
        for change in simulation.algorithm.changes:
            warn('Unsuported algorithm parameter `{}` was ignored.'.format(change.kisao_id), BioSimulatorsWarning)

    # return
    return (
        exec_kisao_id,
        method,
        lambda simulation: alg_props['method_args'](simulation) + method_args,
    )

def exec_simulation(method_name, model, args=None):
    """ Execute a task

    Args:
        method_name (:obj:`str`): name of the :obj:`mpbn` simulation/analysis method
        model (:obj:`py4j.java_gateway.JavaObject`): model
        args (:obj:`list` of :obj:`str`, optional): argument to :obj:`method`

    Returns:
        :obj:`list` of :obj:`dict`: result of :obj:`method` for :obj:`model` and :obj:`args`
    """
    method = getattr(model, method_name)
    if args:
        args_list = [' '.join(args)]
    else:
        args_list = []
    result = method(*args_list)
    return list(result)

def get_variable_target_xpath_ids(variables, model_etree):
    """ Get the SBML-qual id for each XML XPath target of a SED-ML variable

    Args:
        variables (:obj:`list` of :obj:`Variable`): variables of data generators
        model_etree (:obj:`lxml.etree._ElementTree`): element tree for model

    Returns:
        :obj:`dict`: dictionary that maps each variable target to the id of the
            corresponding qualitative species
    """
    namespaces = biosimulators_utils.xml.utils.get_namespaces_for_xml_doc(model_etree)
    sbml_qual_prefix, sbml_qual_uri = biosimulators_utils.model_lang.sbml.utils.get_package_namespace('qual', namespaces)
    return biosimulators_utils.sedml.validation.validate_target_xpaths(
        variables,
        model_etree,
        attr={
            'namespace': {
                'prefix': sbml_qual_prefix,
                'uri': sbml_qual_uri,
            },
            'name': 'id',
        }
    )

def validate_variables(variables, model, model_language, target_xpath_ids, simulation):
    """ Get the result of each SED-ML variable

    Args:
        variables (:obj:`list` of :obj:`Variable`): variables
        model (:obj:`py4j.java_gateway.JavaObject`): bioLQM model
        model_language (:obj:`str`): model language
        target_xpath_ids (:obj:`dict`): dictionary that maps XPaths to the SBML qualitative ids
            of the corresponding objects
        simulation (:obj:`Simulation`): analysis
    """
    component_ids = set(model)
    invalid_variables = []

    for variable in variables:
        if variable.symbol:
            invalid_variables.append('{}: symbol: {}'.format(variable.id, variable.symbol))
        else:
            if model_language == ModelLanguage.SBML.value:
                id = target_xpath_ids[variable.target]
            else:
                id = variable.target

            if id not in component_ids:
                invalid_variables.append('{}: target: {}'.format(variable.id, variable.target))

    if invalid_variables:
        valid_variables = []

        for component_id in component_ids:
            if model_language == ModelLanguage.SBML.value:
                valid_variables.append(
                    "target: /sbml:sbml/sbml:model/qual:listOfQualitativeSpecies/qual:qualitativeSpecies[@id='{}']".format(component_id))
            else:
                valid_variables.append(
                    'target: {}'.format(component_id))

        raise ValueError((
            'The following variables cannot be recorded:\n'
            '  {}\n'
            '\n'
            'Variables with the following symbols and targets can be recorded:\n'
            '  {}'
        ).format(
            '\n  '.join(sorted(invalid_variables)),
            '\n  '.join(sorted(valid_variables)),
        ))

def get_variable_results(variables, model_language, target_xpath_ids, simulation, raw_results):
    """ Get the result of each SED-ML variable

    Args:
        variables (:obj:`list` of :obj:`Variable`): variables
        model_language (:obj:`str`): model language
        target_xpath_ids (:obj:`dict`): dictionary that maps XPaths to the SBML qualitative ids
            of the corresponding objects
        simulation (:obj:`Simulation`): simulation
        raw_results (:obj:`list` of :obj:`dict`): predicted simulatioin states

    Returns:
        :obj:`VariableResults`: result of each SED-ML variable
    """
    n_states = len(raw_results)
    variable_results = VariableResults()
    for variable in variables:
        variable_results[variable.id] = numpy.full((n_states,), numpy.nan)

    for i_state, state in enumerate(raw_results):
        for variable in variables:
            if variable.symbol:
                variable_results[variable.id][i_state] = i_state

            else:
                if model_language == ModelLanguage.SBML.value:
                    id = target_xpath_ids[variable.target]
                else:
                    id = variable.target

                variable_results[variable.id][i_state] = state[id]

    return variable_results

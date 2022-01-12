""" Methods for executing SED tasks and saving their results

:Author: Loïc Paulevé <loic.pauleve@labri.fr>
:Date: 2022-01-12
:Copyright: 2022, Loïc Paulevé
:License: MIT
"""

from .data_model import KISAO_METHOD_MAP
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.config import get_config, Config  # noqa: F401
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog, StandardOutputErrorCapturerLevel  # noqa: F401
from biosimulators_utils.viz.data_model import VizFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults, SedDocumentResults  # noqa: F401
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, ModelAttributeChange,  # noqa: F401
                                                  SteadyStateSimulation, Variable)
from biosimulators_utils.sedml.exec import exec_sed_doc as base_exec_sed_doc
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.utils.core import parse_value, raise_errors_warnings
from kisao.utils import get_preferred_substitute_algorithm_by_ids

from .utils import (
        read_model, exec_simulation, set_up_simulation, get_variable_results,
        get_variable_target_xpath_ids, validate_variables)
from ._version import __version__
import copy
import os
import lxml.etree
import tempfile

__all__ = ['get_simulator_version', 'exec_sedml_docs_in_combine_archive', 'exec_sed_doc', 'exec_sed_task', 'preprocess_sed_task']


def get_simulator_version():
    """ Get the version of MySimulator

    Returns:
        :obj:`str`: version
    """
    return __version__


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            * :obj:`SedDocumentResults`: results
            * :obj:`CombineArchiveLog`: log
    """
    return exec_sedml_docs_in_archive(exec_sed_doc, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      config=config)


def exec_sed_doc(doc, working_dir, base_out_path, rel_out_path=None,
                 apply_xml_model_changes=False,
                 log=None, indent=0, pretty_print_modified_xml_models=False,
                 log_level=StandardOutputErrorCapturerLevel.c, config=None):
    """ Execute the tasks specified in a SED document and generate the specified outputs

    Args:
        doc (:obj:`SedDocument` or :obj:`str`): SED document or a path to SED-ML file which defines a SED document
        working_dir (:obj:`str`): working directory of the SED document (path relative to which models are located)

        base_out_path (:obj:`str`): path to store the outputs

            * CSV: directory in which to save outputs to files
              ``{base_out_path}/{rel_out_path}/{report.id}.csv``
            * HDF5: directory in which to save a single HDF5 file (``{base_out_path}/reports.h5``),
              with reports at keys ``{rel_out_path}/{report.id}`` within the HDF5 file

        rel_out_path (:obj:`str`, optional): path relative to :obj:`base_out_path` to store the outputs
        apply_xml_model_changes (:obj:`bool`, optional): if :obj:`True`, apply any model changes specified in the SED-ML file before
            calling :obj:`task_executer`.
        log (:obj:`SedDocumentLog`, optional): log of the document
        indent (:obj:`int`, optional): degree to indent status messages
        pretty_print_modified_xml_models (:obj:`bool`, optional): if :obj:`True`, pretty print modified XML models
        log_level (:obj:`StandardOutputErrorCapturerLevel`, optional): level at which to log output
        config (:obj:`Config`, optional): BioSimulators common configuration
        simulator_config (:obj:`SimulatorConfig`, optional): tellurium configuration

    Returns:
        :obj:`tuple`:

            * :obj:`ReportResults`: results of each report
            * :obj:`SedDocumentLog`: log of the document
    """
    return base_exec_sed_doc(exec_sed_task, doc, working_dir, base_out_path,
                             rel_out_path=rel_out_path,
                             apply_xml_model_changes=apply_xml_model_changes,
                             log=log,
                             indent=indent,
                             pretty_print_modified_xml_models=pretty_print_modified_xml_models,
                             log_level=log_level,
                             config=config)

def exec_sed_task(task, variables, preprocessed_task=None, log=None, config=None):
    """ Execute a task and save its results

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        preprocessed_task (:obj:`dict`, optional): preprocessed information about the task, including possible
            model changes and variables. This can be used to avoid repeatedly executing the same initialization
            for repeated calls to this method.
        log (:obj:`TaskLog`, optional): log for the task
        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log

    Raises:
        :obj:`NotImplementedError`:

          * Task requires a time course that GINsim doesn't support
          * Task requires an algorithm that GINsim doesn't support
    """
    config = config or get_config()

    if config.LOG and not log:
        log = TaskLog()

    if preprocessed_task is None:
        preprocessed_task = preprocess_sed_task(task, variables, config=config)

    # validate task
    model = task.model
    sim = task.simulation

    # read model
    sim_model = preprocessed_task['model']
    variable_target_sbml_id_map = preprocessed_task['variable_target_sbml_id_map']

    # modify model
    if model.changes:
        if model.language == ModelLanguage.SBML.value:
            allowed_change_types = ()
            # allowed_change_types = (ModelAttributeChange,) # TODO: uncomment when bioLQM recognizes initial state
            raise_errors_warnings(
                validation.validate_model_change_types(model.changes, allowed_change_types),
                error_summary='Changes for model `{}` are not supported.'.format(model.id))

            model_etree = preprocessed_task['model_etree']

            model = copy.deepcopy(model)
            for change in model.changes:
                change.new_value = str(int(float(change.new_value)))

            apply_changes_to_xml_model(model, model_etree, sed_doc=None, working_dir=None)

            model_file, model_filename = tempfile.mkstemp(suffix='.xml')
            os.close(model_file)

            model_etree.write(model_filename,
                              xml_declaration=True,
                              encoding="utf-8",
                              standalone=False,
                              pretty_print=False)

            sim_model = read_model(model_filename, model.language)

            os.remove(model_filename)

        else:
            raise_errors_warnings(
                validation.validate_model_change_types(model.changes, ()),
                error_summary='Changes for model `{}` are not supported.'.format(model.id))

    alg_method_name = preprocessed_task['algorithm_method_name']
    alg_method_args = preprocessed_task['algorithm_method_args'](sim)

    # run simulation
    raw_results = exec_simulation(alg_method_name, sim_model, alg_method_args)

    # transform results
    variable_results = get_variable_results(variables, model.language, variable_target_sbml_id_map, sim, raw_results)

    # log action
    if config.LOG:
        log.algorithm = preprocessed_task['algorithm_kisao_id']
        log.simulator_details = {
            'method': alg_method_name,
            'arguments': alg_method_args,
        }

    ############################
    # return the result of each variable and log
    return variable_results, log

def preprocess_sed_task(task, variables, config=None):
    """ Preprocess a SED task, including its possible model changes and variables. This is useful for avoiding
    repeatedly initializing tasks on repeated calls of :obj:`exec_sed_task`.

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`object`: preprocessed information about the task
    """
    config = config or get_config()

    # validate task
    model = task.model
    sim = task.simulation

    if config.VALIDATE_SEDML:
        raise_errors_warnings(
            validation.validate_task(task),
            error_summary='Task `{}` is invalid.'.format(task.id))
        raise_errors_warnings(
            validation.validate_model_language(
                model.language,
                (ModelLanguage.SBML, ModelLanguage.ZGINML)),
            error_summary='Language for model `{}` is not supported.'.format(model.id))

        if model.language == ModelLanguage.SBML.value:
            allowed_change_types = (ModelAttributeChange,)
        else:
            allowed_change_types = ()
        raise_errors_warnings(
            validation.validate_model_change_types(model.changes, allowed_change_types),
            error_summary='Changes for model `{}` are not supported.'.format(model.id))
        raise_errors_warnings(
            *validation.validate_model_changes(model),
            error_summary='Changes for model `{}` are invalid.'.format(model.id))
        raise_errors_warnings(
            validation.validate_simulation_type(sim, (SteadyStateSimulation,)),
            error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))

    if not os.path.isfile(model.source):
        raise FileNotFoundError('`{}` is not a file.'.format(model.source))

    if model.language == ModelLanguage.SBML.value:
        model_etree = lxml.etree.parse(model.source)
        variable_target_sbml_id_map = get_variable_target_xpath_ids(variables, model_etree)
    else:
        model_etree = None
        variable_target_sbml_id_map = None

    # validate model
    if config.VALIDATE_SEDML_MODELS:
        raise_errors_warnings(*validation.validate_model(model, [], working_dir='.'),
                              error_summary='Model `{}` is invalid.'.format(model.id),
                              warning_summary='Model `{}` may be invalid.'.format(model.id))

    # read model
    sim_model = read_model(model.source, model.language)

    # validate variables
    validate_variables(variables, sim_model, model.language, variable_target_sbml_id_map, sim)

    # setup simulation
    alg_kisao_id, method_name, method_args = set_up_simulation(sim, config=config)

    ############################
    # return preprocessed information
    return {
        'model': sim_model,
        'model_etree': model_etree,
        'variable_target_sbml_id_map': variable_target_sbml_id_map,
        'algorithm_kisao_id': alg_kisao_id,
        'algorithm_method_name': method_name,
        'algorithm_method_args': method_args,
    }

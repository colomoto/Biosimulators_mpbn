""" Tests of the command-line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-07-07
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_mpbn import __main__
from biosimulators_mpbn import core
from biosimulators_mpbn.data_model import KISAO_METHOD_MAP
from biosimulators_mpbn.utils import read_model
from biosimulators_utils.combine import data_model as combine_data_model
from biosimulators_utils.combine.io import CombineArchiveWriter
from biosimulators_utils.config import get_config
from biosimulators_utils.report import data_model as report_data_model
from biosimulators_utils.report.io import ReportReader
from biosimulators_utils.simulator.exec import exec_sedml_docs_in_archive_with_containerized_simulator
from biosimulators_utils.simulator.specs import gen_algorithms_from_specs
from biosimulators_utils.sedml import data_model as sedml_data_model
from biosimulators_utils.sedml.io import SedmlSimulationWriter
from biosimulators_utils.sedml.utils import append_all_nested_children_to_doc
from biosimulators_utils.warnings import BioSimulatorsWarning
from kisao.exceptions import AlgorithmCannotBeSubstitutedException
from unittest import mock
import datetime
import dateutil.tz
import json
import numpy
import numpy.testing
import os
import random
import shutil
import tempfile
import unittest
import yaml


class CliTestCase(unittest.TestCase):
    EXAMPLE_SBML_MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'example-model.xml')
    EXAMPLE_ZGINML_MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'SuppMat_Model_Master_Model.zginml')
    SPECIFICATIONS_FILENAME = os.path.join(os.path.dirname(__file__), '..', 'biosimulators.json')
    DOCKER_IMAGE = 'ghcr.io/biosimulators/biosimulators_mpbn/mpbn:latest'
    NAMESPACES = {
        None: 'http://sed-ml.org/sed-ml/level1/version3',
        'sbml': 'http://www.sbml.org/sbml/level3/version1/core',
        'qual': 'http://www.sbml.org/sbml/level3/version1/qual/version1',
    }

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_exec_sbml_sed_task_successfully(self):
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_SBML_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.SBML.value,
            ),
            simulation=sedml_data_model.SteadyStateSimulation(
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000660',
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='G0',
                target="/sbml:sbml/sbml:model/qual:listOfQualitativeSpecies/qual:qualitativeSpecies[@qual:id='G0']/@level",
                target_namespaces=self.NAMESPACES,
                task=task),
            sedml_data_model.Variable(
                id='G1',
                target="/sbml:sbml/sbml:model/qual:listOfQualitativeSpecies/qual:qualitativeSpecies[@qual:id='G1']/@level",
                target_namespaces=self.NAMESPACES,
                task=task),
        ]

        task.simulation.algorithm.kisao_id = 'KISAO_0000660'
        variable_results, log = core.exec_sed_task(task, variables)
        self.assertEqual(set(variable_results.keys()), set(['G0', 'G1']))
        for variable_result in variable_results.values():
            self.assertFalse(numpy.any(numpy.isnan(variable_result)))

        # check that log can be serialized to JSON
        json.dumps(log.to_json())

        log.out_dir = self.dirname
        log.export()
        with open(os.path.join(self.dirname, get_config().LOG_PATH), 'rb') as file:
            log_data = yaml.load(file, Loader=yaml.Loader)
        json.dumps(log_data)


    def test_exec_zginml_sed_task_successfully(self):
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=self.EXAMPLE_ZGINML_MODEL_FILENAME,
                language=sedml_data_model.ModelLanguage.ZGINML.value,
            ),
            simulation=sedml_data_model.SteadyStateSimulation(
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000661',
                ),
            ),
        )

        variables = [
            sedml_data_model.Variable(
                id='AKT1',
                target="AKT1",
                task=task),
            sedml_data_model.Variable(
                id='DKK1',
                target="DKK1",
                task=task),
        ]

        # synchronous method
        task.simulation.algorithm.kisao_id = 'KISAO_0000661'
        variable_results, log = core.exec_sed_task(task, variables)
        self.assertEqual(set(variable_results.keys()), set(['AKT1', 'DKK1']))
        for variable_result in variable_results.values():
            self.assertFalse(numpy.any(numpy.isnan(variable_result)))

        # check that log can be serialized to JSON
        json.dumps(log.to_json())

        log.out_dir = self.dirname
        log.export()
        with open(os.path.join(self.dirname, get_config().LOG_PATH), 'rb') as file:
            log_data = yaml.load(file, Loader=yaml.Loader)
        json.dumps(log_data)

        # error handling
        preprocessed_task = core.preprocess_sed_task(task, variables)
        task.model.changes.append(None)
        with self.assertRaisesRegex(ValueError, 'are not supported'):
            core.exec_sed_task(task, variables, preprocessed_task=preprocessed_task)

    def test_exec_sedml_docs_in_combine_archive_successfully(self):
        doc, archive_filename = self._build_combine_archive()

        out_dir = os.path.join(self.dirname, 'out')

        config = get_config()
        config.REPORT_FORMATS = [report_data_model.ReportFormat.h5]
        config.BUNDLE_OUTPUTS = True
        config.KEEP_INDIVIDUAL_OUTPUTS = True

        _, log = core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=config)
        if log.exception:
            raise log.exception

        self._assert_combine_archive_outputs(doc, out_dir)

    def _build_combine_archive(self, algorithm=None):
        doc = self._build_sed_doc(algorithm=algorithm)

        archive_dirname = os.path.join(self.dirname, 'archive')
        if not os.path.isdir(archive_dirname):
            os.mkdir(archive_dirname)

        model_filename = os.path.join(archive_dirname, 'model.xml')
        shutil.copyfile(self.EXAMPLE_SBML_MODEL_FILENAME, model_filename)

        sim_filename = os.path.join(archive_dirname, 'sim.sedml')
        SedmlSimulationWriter().run(doc, sim_filename)

        archive = combine_data_model.CombineArchive(
            contents=[
                combine_data_model.CombineArchiveContent(
                    'model.xml', combine_data_model.CombineArchiveContentFormat.SBML.value),
                combine_data_model.CombineArchiveContent(
                    'sim.sedml', combine_data_model.CombineArchiveContentFormat.SED_ML.value),
            ],
        )
        archive_filename = os.path.join(self.dirname, 'archive.omex')

        CombineArchiveWriter().run(archive, archive_dirname, archive_filename)

        shutil.copy(archive_filename, f"/tmp/archive_{random.randint(1000,2000)}.omex")

        return (doc, archive_filename)

    def _build_sed_doc(self, algorithm=None):
        if algorithm is None:
            algorithm = sedml_data_model.Algorithm(
                kisao_id='KISAO_0000660',
            )

        doc = sedml_data_model.SedDocument()
        doc.models.append(sedml_data_model.Model(
            id='model',
            source='model.xml',
            language=sedml_data_model.ModelLanguage.SBML.value,
        ))
        doc.simulations.append(sedml_data_model.SteadyStateSimulation(
            id='sim_steady_state',
            algorithm=algorithm,
        ))

        doc.tasks.append(sedml_data_model.Task(
            id='task_1',
            model=doc.models[0],
            simulation=doc.simulations[0],
        ))

        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_G1',
            variables=[
                sedml_data_model.Variable(
                    id='var_G1',
                    target="/sbml:sbml/sbml:model/qual:listOfQualitativeSpecies/qual:qualitativeSpecies[@qual:id='G1']",
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_G1',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_G2',
            variables=[
                sedml_data_model.Variable(
                    id='var_G2',
                    target="/sbml:sbml/sbml:model/qual:listOfQualitativeSpecies/qual:qualitativeSpecies[@qual:id='G2']/@level",
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_G2',
        ))

        doc.outputs.append(sedml_data_model.Report(
            id='report',
            data_sets=[
                sedml_data_model.DataSet(id='data_set_G1', label='G1', data_generator=doc.data_generators[0]),
                sedml_data_model.DataSet(id='data_set_G2', label='G2', data_generator=doc.data_generators[1]),
            ],
        ))
        if isinstance(doc.simulations[0], sedml_data_model.UniformTimeCourseSimulation):
            doc.data_generators.insert(
                0,
                sedml_data_model.DataGenerator(
                    id='data_gen_time',
                    variables=[
                        sedml_data_model.Variable(
                            id='var_time',
                            symbol=sedml_data_model.Symbol.time,
                            task=doc.tasks[0],
                        ),
                    ],
                    math='var_time',
                ))
            doc.outputs[0].data_sets.insert(
                0,
                sedml_data_model.DataSet(id='data_set_time', label='Time', data_generator=doc.data_generators[0]))

        append_all_nested_children_to_doc(doc)

        return doc

    def _assert_combine_archive_outputs(self, doc, out_dir):
        self.assertEqual(set(['reports.h5']).difference(set(os.listdir(out_dir))), set())

        report = ReportReader().run(doc.outputs[0], out_dir, 'sim.sedml/report', format=report_data_model.ReportFormat.h5)

        self.assertEqual(sorted(report.keys()), sorted([d.id for d in doc.outputs[0].data_sets]))

        sim = doc.tasks[0].simulation
        if isinstance(sim, sedml_data_model.UniformTimeCourseSimulation):
            self.assertEqual(len(report[doc.outputs[0].data_sets[0].id]), sim.number_of_points + 1)

        for data_set_result in report.values():
            self.assertFalse(numpy.any(numpy.isnan(data_set_result)))

        if isinstance(sim, sedml_data_model.UniformTimeCourseSimulation):
            self.assertIn('data_set_time', report)
            numpy.testing.assert_allclose(report[doc.outputs[0].data_sets[0].id],
                                          numpy.linspace(sim.output_start_time, sim.output_end_time, sim.number_of_points + 1))
        else:
            self.assertNotIn('data_set_time', report)

    def test_exec_sedml_docs_in_combine_archive_with_all_algorithms(self):
        for alg in gen_algorithms_from_specs(self.SPECIFICATIONS_FILENAME).values():
            doc, archive_filename = self._build_combine_archive(algorithm=alg)
            out_dir = os.path.join(self.dirname, alg.kisao_id)

            config = get_config()
            config.REPORT_FORMATS = [report_data_model.ReportFormat.h5]
            config.BUNDLE_OUTPUTS = True
            config.KEEP_INDIVIDUAL_OUTPUTS = True

            _, log = core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=config)
            if log.exception:
                raise log.exception

            self._assert_combine_archive_outputs(doc, out_dir)

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: ')

    def test_exec_sedml_docs_in_combine_archive_with_cli(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        env = self._get_combine_archive_exec_env()

        with mock.patch.dict(os.environ, env):
            with __main__.App(argv=['-i', archive_filename, '-o', out_dir]) as app:
                app.run()

        self._assert_combine_archive_outputs(doc, out_dir)

    def _get_combine_archive_exec_env(self):
        return {
            'REPORT_FORMATS': 'h5'
        }

    def test_exec_sedml_docs_in_combine_archive_with_docker_image(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        docker_image = self.DOCKER_IMAGE
        env = self._get_combine_archive_exec_env()

        exec_sedml_docs_in_archive_with_containerized_simulator(
            archive_filename, out_dir, docker_image, environment=env, pull_docker_image=False)

        self._assert_combine_archive_outputs(doc, out_dir)

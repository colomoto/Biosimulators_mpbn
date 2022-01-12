from biosimulators_mpbn.utils import (get_variable_target_xpath_ids,
                                        validate_variables, read_model, set_up_simulation,
                                        exec_simulation, get_variable_results)
from biosimulators_utils.sedml.data_model import (ModelLanguage, SteadyStateSimulation,
                                                  Algorithm, AlgorithmParameterChange,
                                                  Variable, Symbol)
from biosimulators_utils.warnings import BioSimulatorsWarning
from kisao.exceptions import AlgorithmCannotBeSubstitutedException
from kisao.warnings import AlgorithmSubstitutedWarning
from unittest import mock
import collections
import lxml.etree
import numpy.testing
import os
import unittest


class UtilsTestCase(unittest.TestCase):
    def test_get_variable_target_xpath_ids(self):
        with mock.patch('lxml.etree.parse', return_value=None):
            with mock.patch('biosimulators_utils.xml.utils.get_namespaces_for_xml_doc', return_value={'qual': None}):
                with mock.patch('biosimulators_utils.model_lang.sbml.utils.get_package_namespace', return_value=('qual', None)):
                    with mock.patch('biosimulators_utils.sedml.validation.validate_target_xpaths', return_value={'x': 'X'}):
                        self.assertEqual(get_variable_target_xpath_ids([Variable(target='x')], None), {'x': 'X'})

    def test_read_model(self):
        filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'SuppMat_Model_Master_Model.zginml')
        read_model(filename, ModelLanguage.ZGINML)

        filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'example-model.sbml')
        read_model(filename, ModelLanguage.SBML)

        filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'example-model.xml')
        read_model(filename, ModelLanguage.SBML)

        with self.assertRaises(FileNotFoundError):
            read_model('not a file', ModelLanguage.SBML)

        filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'invalid.zginml')
        with self.assertRaises(ValueError):
            read_model(filename, ModelLanguage.ZGINML)

        filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'regulatoryGraph.ginml')
        with self.assertRaises(ValueError):
            read_model(filename, ModelLanguage.GINML)

    def test_set_up_simulation(self):
        simulation = SteadyStateSimulation(algorithm=Algorithm())
        simulation.algorithm.kisao_id = 'KISAO_0000660'
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            kisao_id, method, method_args = set_up_simulation(simulation)
        self.assertEqual(kisao_id, 'KISAO_0000660')
        self.assertEqual(method, 'fixedpoints')
        self.assertEqual(method_args(simulation), [])

        simulation.algorithm.kisao_id = 'KISAO_0000661'
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            kisao_id, method, method_args = set_up_simulation(simulation)
        self.assertEqual(kisao_id, 'KISAO_0000661')
        self.assertEqual(method, 'attractors')
        self.assertEqual(method_args(simulation), [])

    def test_exec_simulation(self):
        filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'SuppMat_Model_Master_Model.zginml')
        model = read_model(filename, ModelLanguage.ZGINML)
        raw_results = exec_simulation('fixedpoints', model)
        self.assertEqual(len(raw_results), 9)

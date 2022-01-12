from biosimulators_mpbn.data_model import KISAO_METHOD_MAP
import json
import os
import unittest


class DataModelTestCase(unittest.TestCase):
    def test_data_model_matches_specs(self):
        specs_filename = os.path.join(os.path.dirname(__file__), '..', 'biosimulators.json')
        with open(specs_filename, 'r') as file:
            specs = json.load(file)

        self.assertEqual(set(KISAO_METHOD_MAP.keys()),
                         set(alg_specs['kisaoId']['id'] for alg_specs in specs['algorithms']))

        for alg_specs in specs['algorithms']:
            alg_props = KISAO_METHOD_MAP[alg_specs['kisaoId']['id']]
            self.assertEqual(set(alg_props['parameters'].keys()),
                             set(param_specs['kisaoId']['id'] for param_specs in alg_specs['parameters']))

            for param_specs in alg_specs['parameters']:
                param_props = alg_props['parameters'][param_specs['kisaoId']['id']]
                self.assertEqual(param_props['type'].value, param_specs['type'])

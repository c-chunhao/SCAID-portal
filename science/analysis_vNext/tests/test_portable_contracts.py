"""Synthetic tests: no real assay, cell or individual records are required."""
from pathlib import Path
import copy
import importlib.util
import json
import os
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gate = module('science_gate', 'release_scientific_gate.py')
paths = module('science_paths', 'runtime_paths.py')
composition = module('coverage_composition', 'statistics/recompute_atlas_composition.py')


class PortableContracts(unittest.TestCase):
    def setUp(self):
        self.record = json.loads((ROOT / 'schemas/release_evidence.example.json').read_text())
        for key, value in list(self.record.items()):
            if isinstance(value, bool):
                self.record[key] = True
        self.record['same_donor_overlap_with_validation'] = False
        self.record['source_input_sha256'] = 'a' * 64

    def test_composition_identity_delimiter_collision_is_rejected(self):
        import pandas as pd
        frame = pd.DataFrame({'OldCells':['a|b'], 'Sample':['c'], 'Dataset':['d']})
        with self.assertRaises(ValueError):
            composition.identity_key(frame)
        frame = pd.DataFrame({'OldCells':['a'], 'Sample':['b|c'], 'Dataset':['d']})
        with self.assertRaises(ValueError):
            composition.identity_key(frame)
        clean = pd.DataFrame({'OldCells':['a'], 'Sample':['b'], 'Dataset':['d']})
        self.assertEqual(composition.identity_key(clean).iloc[0], 'a|b|d')

    def test_complete_boolean_evidence_permits_completeness_only(self):
        result = gate.evaluate(self.record)
        self.assertTrue(result['new_release_eligible'])
        self.assertTrue(result['donor_disease_contrast_eligible'])
        self.assertTrue(result['evidence_completeness_only'])

    def test_false_string_does_not_pass_release(self):
        self.record['counts_schema_verified'] = 'false'
        self.assertFalse(gate.evaluate(self.record)['new_release_eligible'])

    def test_true_string_does_not_pass_release(self):
        self.record['annotation_review_completed'] = 'true'
        self.assertFalse(gate.evaluate(self.record)['new_release_eligible'])

    def test_numeric_one_does_not_pass_release(self):
        self.record['capture_mapping_verified'] = 1
        self.assertFalse(gate.evaluate(self.record)['new_release_eligible'])

    def test_nonhex_64_characters_are_not_SHA256(self):
        self.record['source_input_sha256'] = 'z' * 64
        self.assertFalse(gate.evaluate(self.record)['new_release_eligible'])

    def test_donor_field_is_not_promoted_from_truthy_string(self):
        self.record['independent_individuals'] = 'false'
        result = gate.evaluate(self.record)
        self.assertTrue(result['new_release_eligible'])
        self.assertFalse(result['donor_disease_contrast_eligible'])

    def test_validation_overlap_fails_donor_inference(self):
        self.record['same_donor_overlap_with_validation'] = True
        self.assertFalse(gate.evaluate(self.record)['donor_disease_contrast_eligible'])

    def test_relative_revision_resolved_exactly(self):
        with tempfile.TemporaryDirectory() as tmp:
            registry = Path(tmp) / 'assets.json'
            rel = 'revisions/2026-10-05/RNA/example.h5'
            registry.write_text(json.dumps({'datasets': [{'h5_paths': {'RNA': rel}}]}))
            env = {'SCAID_ASSET_REGISTRY': str(registry), 'SCAID_DATA_H5_ROOT': str(Path(tmp) / 'assays')}
            with patch.dict(os.environ, env, clear=True):
                result = paths.registered_datasets()
            self.assertEqual(result[0]['h5_paths']['RNA'], str(Path(tmp) / 'assays' / rel))

    def test_relative_registry_requires_explicit_local_base(self):
        with tempfile.TemporaryDirectory() as tmp:
            registry = Path(tmp) / 'assets.json'
            registry.write_text(json.dumps({'datasets': [{'h5_paths': {'RNA': 'RNA/example.h5'}}]}))
            with patch.dict(os.environ, {'SCAID_ASSET_REGISTRY': str(registry)}, clear=True):
                with self.assertRaises(ValueError):
                    paths.registered_datasets()

    def test_parent_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            registry = Path(tmp) / 'assets.json'
            registry.write_text(json.dumps({'datasets': [{'h5_paths': {'RNA': '../example.h5'}}]}))
            with patch.dict(os.environ, {'SCAID_ASSET_REGISTRY': str(registry), 'SCAID_DATA_H5_ROOT': tmp}, clear=True):
                with self.assertRaises(ValueError):
                    paths.registered_datasets()


if __name__ == '__main__':
    unittest.main()

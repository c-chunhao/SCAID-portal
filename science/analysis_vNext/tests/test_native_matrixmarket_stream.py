"""Compile native source temporarily and test eight anonymous stream fixtures."""
from pathlib import Path
import gzip
import json
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('g++'), 'g++ is required for native scanner tests')
class NativeMatrixMarket(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.tmp.name)
        cls.executable = cls.directory / 'audit_matrixmarket_gz'
        subprocess.run(['g++', '-O3', '-std=c++17',
                        str(ROOT / 'new_data/audit_matrixmarket_gz.cpp'),
                        '-lz', '-o', str(cls.executable)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_complete_stream_scope_and_failures(self):
        cases = {
            'valid': '%%MatrixMarket matrix coordinate integer general\n2 3 2\n1 1 4\n2 3 7\n',
            'fractional': '%%MatrixMarket matrix coordinate real general\n2 3 2\n1 1 4.5\n2 3 7\n',
            'missing_record': '%%MatrixMarket matrix coordinate integer general\n2 3 2\n1 1 4\n',
            'out_of_bounds': '%%MatrixMarket matrix coordinate integer general\n2 3 1\n3 1 4\n',
        }
        for case, source in cases.items():
            for compressed in [False, True]:
                with self.subTest(case=case, gzip=compressed):
                    path = self.directory / (case + ('.mtx.gz' if compressed else '.mtx'))
                    if compressed:
                        with gzip.open(path, 'wt') as handle:
                            handle.write(source)
                    else:
                        path.write_text(source)
                    run = subprocess.run([str(self.executable), str(path)],
                                         capture_output=True, text=True)
                    if case in ['valid', 'fractional']:
                        self.assertEqual(run.returncode, 0)
                        report = json.loads(run.stdout)
                        self.assertTrue(report['schema_pass'])
                        self.assertEqual(report['numeric_entries_scanned'], 2)
                        self.assertEqual(report['full_numeric'],
                                         'integer_nonnegative' if case == 'valid' else 'noninteger_or_invalid')
                    else:
                        self.assertNotEqual(run.returncode, 0)


if __name__ == '__main__':
    unittest.main()

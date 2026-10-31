"""Source provenance and reproducible archive contracts."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('gene_build',ROOT/'scripts/build_release.py')
build=importlib.util.module_from_spec(spec); spec.loader.exec_module(build)


class PackagingTests(unittest.TestCase):
    def test_source_provenance(self):
        record=json.loads((ROOT/'provenance/source_snapshot.json').read_text())
        self.assertIsNone(record['upstream_git_commit'])
        self.assertEqual(len({x['path'] for x in record['files']}),len(record['files']))
        for item in record['files']:
            p=ROOT/item['path']
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),item['packaged_sha256'],item['path'])
            self.assertEqual(p.stat().st_size,item['packaged_bytes'])
    def test_archive_reproducibility(self):
        with tempfile.TemporaryDirectory(prefix='gene_release_') as directory:
            a,b=Path(directory)/'a',Path(directory)/'b'; build.build(a); build.build(b)
            for p in a.iterdir(): self.assertEqual(p.read_bytes(),(b/p.name).read_bytes())
            with zipfile.ZipFile(next(a.glob('*.zip'))) as archive:
                names=archive.namelist()
                self.assertTrue(any(n.endswith('/MANIFEST.sha256') for n in names))
                self.assertFalse(any('__pycache__' in n or n.endswith('.sif') for n in names))
    def test_disabled_recovery_and_native_hints(self):
        import configparser
        cfg=configparser.ConfigParser(); cfg.read(ROOT/'config/braker.parameters.ini')
        self.assertEqual(cfg.getint('PARAMS','run_best_by_compleasm'),0)
        self.assertEqual(cfg.getint('PARAMS','use_compleasm_hints'),1)
        self.assertEqual(cfg.getint('PARAMS','skip_busco'),1)


if __name__=='__main__': unittest.main()

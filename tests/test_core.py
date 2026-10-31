"""Portable unit tests: no prediction binaries, data downloads or containers."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'pipeline'))
import pipeline as core
import filter_engine as engine


def row(label, c, genes, proteins=None, intron=.8):
    return dict(label=label, annotation=dict(genes=genes, proteins=proteins or genes),
                parameters=dict(intron_support=intron),
                assessments=dict(arachnida=dict(N=100, C=c)))


class SelectionTests(unittest.TestCase):
    def test_exact_panel(self):
        self.assertEqual(core.PANEL, [('baseline', 1.0, 2, True), ('A', .9, 2, True),
            ('B', .8, 2, True), ('C', .8, 1, True), ('D', .8, 1, False),
            ('E', .75, 2, True), ('F', .7, 2, True), ('G', .6, 2, True), ('H', .5, 2, True)])

    def test_reached_target_prefers_compact_not_highest_c(self):
        selected, _ = core.choose([row('A', 95, 30), row('B', 99, 40)], dict(N=100, C=95))
        self.assertEqual(selected['label'], 'A')

    def test_unmet_target_prefers_c(self):
        selected, reason = core.choose([row('A', 90, 30), row('B', 94, 40)], dict(N=100, C=95))
        self.assertEqual(selected['label'], 'B'); self.assertIn('No setting reached', reason)

    def test_ties(self):
        rows = [row('A', 95, 30, 40), row('B', 95, 30, 35, .7), row('C', 95, 30, 35, .8)]
        self.assertEqual(core.choose(rows, dict(N=100, C=95))[0]['label'], 'C')
        self.assertEqual(core.choose([row('B', 95, 30), row('A', 95, 30)], dict(N=100, C=95))[0]['label'], 'A')

    def test_lineage_size_mismatch(self):
        with self.assertRaises(AssertionError): core.choose([row('A', 95, 30)], dict(N=99, C=95))

    def test_mode_candidate_sets(self):
        engine.set_mode('ep'); self.assertIn('GeneMark-EP/genemark.gtf', engine.INPUTS)
        self.assertIn('hintsfile_iter2.gff', engine.INPUTS)
        engine.set_mode('etp'); self.assertIn('GeneMark-ETP/training.gtf', engine.INPUTS)
        self.assertNotIn('hintsfile_iter2.gff', engine.INPUTS)


@unittest.skipUnless(sys.platform.startswith('linux'), 'POSIX ownership and cleanup checks')
class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='gene_test_')
        self.root = Path(self.tmp.name).resolve()
        (self.root / 'a.txt').write_text('a\n')
    def tearDown(self): self.tmp.cleanup()

    def test_checksum_tamper(self):
        engine.write_manifest(self.root); self.assertEqual(engine.verify_manifest(self.root), 1)
        (self.root / 'a.txt').write_text('changed')
        with self.assertRaises(AssertionError): engine.verify_manifest(self.root)

    def test_unlisted_files_fail(self):
        engine.write_manifest(self.root); (self.root / 'b.txt').write_text('b')
        with self.assertRaises(AssertionError): engine.verify_manifest(self.root)

    def test_symlinks_fail(self):
        (self.root / 'b.txt').symlink_to(self.root / 'a.txt')
        with self.assertRaises(AssertionError): engine.write_manifest(self.root)

    def test_checksum_path_escape_fails(self):
        (self.root / 'SHA256SUMS').write_text('0'*64 + '  ../outside\n')
        with self.assertRaises(AssertionError): engine.verify_manifest(self.root)

    def test_owner_guard(self):
        work = self.root / 'work'; work.mkdir(); (work / 'unrelated').write_text('keep')
        with self.assertRaises(AssertionError): core.own_temporary(work, self.root / 'out')
        self.assertTrue((work / 'unrelated').exists())


if __name__ == '__main__': unittest.main()

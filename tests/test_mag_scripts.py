#!/usr/bin/env python3
"""Tests for MAG-module helper scripts."""

import csv
import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = str(REPO / "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


def load_script(name):
    path = REPO / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# extract_genomic_context
# ---------------------------------------------------------------------------

GFF_CONTENT = """\
##gff-version 3
##sequence-region contig1 1 200000
contig1\tProdigal\tCDS\t1000\t2000\t.\t+\t.\tID=PROKKA_00001;gene=intA;product=integrase
contig1\tProdigal\tCDS\t3000\t4000\t.\t+\t.\tID=PROKKA_00002;gene=clbA_prot;product=colibactin protein A
contig1\tProdigal\tCDS\t100000\t101000\t.\t-\t.\tID=PROKKA_00003;gene=tnpA;product=transposase
contig2\tProdigal\tCDS\t500\t1500\t.\t+\t.\tID=PROKKA_00004;gene=gyrB;product=gyrase subunit B
"""

GFF_WITH_GENE_FEATURES = """\
##gff-version 3
##sequence-region contig1 1 200000
contig1\tProkka\tgene\t1000\t2000\t.\t+\t.\tID=gene_1
contig1\tProkka\tCDS\t1000\t2000\t.\t+\t.\tID=PROKKA_00001;gene=intA;product=integrase
contig1\tProkka\tgene\t3000\t4000\t.\t+\t.\tID=gene_2
contig1\tProkka\tCDS\t3000\t4000\t.\t+\t.\tID=PROKKA_00002;gene=clbA;product=colibactin A
"""

TBLOUT_CONTENT = """\
#                                                               --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias  exp  dom  seq  from    to  from    to  from    to  acc description of target
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ ----- ---- ---- ---- ----- ----- ----- ----- ----- ----- ---- ---------------------
contig1_2            -          clbA                 -            1.2e-10  35.0   0.0   1.5e-10   34.8   0.0   1.0    1    1     1   200    10   205    8   207  0.95 colibactin A
"""


class TestExtractGenomicContext(unittest.TestCase):
    def setUp(self):
        self.mod = load_script("extract_genomic_context")

    def _write(self, tmp, name, content):
        p = os.path.join(tmp, name)
        with open(p, "w") as f:
            f.write(content)
        return p

    def test_parse_prokka_gff_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = self._write(tmp, "test.gff", GFF_CONTENT)
            genes = self.mod.parse_prokka_gff(gff)
        self.assertEqual(len(genes), 4)

    def test_parse_prokka_gff_excludes_gene_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = self._write(tmp, "test.gff", GFF_WITH_GENE_FEATURES)
            genes = self.mod.parse_prokka_gff(gff)
        self.assertEqual(len(genes), 2)

    def test_parse_prokka_gff_gene_names_from_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = self._write(tmp, "test.gff", GFF_WITH_GENE_FEATURES)
            genes = self.mod.parse_prokka_gff(gff)
        self.assertEqual(genes[0]["locus_tag"], "PROKKA_00001")
        self.assertEqual(genes[0]["gene"], "intA")

    def test_parse_hmmsearch_tblout(self):
        with tempfile.TemporaryDirectory() as tmp:
            tbl = self._write(tmp, "test.tblout", TBLOUT_CONTENT)
            hits = self.mod.parse_hmmsearch_tblout(tbl, evalue_threshold=1e-5)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["clb_gene"], "clbA")
        self.assertEqual(hits[0]["contig"], "contig1")

    def test_is_integrase(self):
        self.assertTrue(self.mod.is_integrase({"gene": "intA", "product": "integrase"}))
        self.assertFalse(self.mod.is_integrase({"gene": "gyrB", "product": "gyrase"}))

    def test_is_transposase(self):
        self.assertTrue(self.mod.is_transposase({"gene": "tnpA", "product": "transposase"}))
        self.assertFalse(self.mod.is_transposase({"gene": "gyrB", "product": "gyrase"}))

    def test_get_flanking_within_window(self):
        genes = [
            {"contig": "c1", "start": 1000, "end": 2000, "locus_tag": "L1", "gene": "g1", "product": ""},
            {"contig": "c1", "start": 80000, "end": 81000, "locus_tag": "L2", "gene": "g2", "product": ""},
            {"contig": "c2", "start": 1000, "end": 2000, "locus_tag": "L3", "gene": "g3", "product": ""},
        ]
        result = self.mod.get_flanking(genes, "c1", 5000, window=50000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["locus_tag"], "L1")

    def test_end_to_end_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = self._write(tmp, "bin.gff", GFF_CONTENT)
            tbl = self._write(tmp, "bin.tblout", TBLOUT_CONTENT)
            out = os.path.join(tmp, "out.tsv")
            import sys
            sys.argv = ["extract_genomic_context.py",
                        "--gff", gff, "--tblout", tbl,
                        "--evalue", "1e-5", "--out", out]
            self.mod.main()
            rows = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["clb_gene"], "clbA")
        # integrase (PROKKA_00001) is within 50kb of contig1_2
        self.assertEqual(rows[0]["has_integrase"], "True")
        self.assertEqual(rows[0]["has_transposase"], "False")


# ---------------------------------------------------------------------------
# build_mag_summary
# ---------------------------------------------------------------------------

CHECKM2_CONTENT = """\
Name\tCompleteness\tContamination\tGenome_size
bin_001\t92.5\t1.2\t4500000
bin_002\t45.0\t8.0\t2100000
"""

GTDBTK_CONTENT = """\
user_genome\tclassification
bin_001\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
bin_002\td__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus acidophilus
"""

TBLOUT_PKS_CONTENT = """\
#
bin_001_1\t-\tclbA\t-\t1e-10\t35.0\t0.0\t1.5e-10\t34.8\t0.0\t1.0\t1\t1\t1\t200\t10\t205\t8\t207\t0.95\thit
bin_001_2\t-\tclbB\t-\t2e-08\t30.0\t0.0\t2.5e-08\t29.8\t0.0\t1.0\t1\t1\t1\t200\t10\t205\t8\t207\t0.95\thit
"""

CONTEXT_CONTENT = """\
locus_tag\tclb_gene\tevalue\tcontig\thas_integrase\thas_transposase\tflanking_genes
bin_001_1\tclbA\t1e-10\tbin_001\tTrue\tFalse\tgyrB;fliA
"""


class TestBuildMagSummary(unittest.TestCase):
    def setUp(self):
        self.mod = load_script("build_mag_summary")

    def _write(self, tmp, name, content):
        p = os.path.join(tmp, name)
        with open(p, "w") as f:
            f.write(content)
        return p

    def test_parse_checkm2(self):
        with tempfile.TemporaryDirectory() as tmp:
            tsv = self._write(tmp, "quality.tsv", CHECKM2_CONTENT)
            result = self.mod.parse_checkm2(tsv)
        self.assertAlmostEqual(result["bin_001"]["completeness"], 92.5)
        self.assertAlmostEqual(result["bin_001"]["contamination"], 1.2)

    def test_parse_gtdbtk(self):
        with tempfile.TemporaryDirectory() as tmp:
            tsv = self._write(tmp, "gtdbtk.tsv", GTDBTK_CONTENT)
            result = self.mod.parse_gtdbtk(tsv)
        self.assertIn("Enterobacterales", result["bin_001"])
        self.assertNotIn("Enterobacterales", result["bin_002"])

    def test_unexpected_taxon_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            gtdbtk = self._write(tmp, "gtdbtk.tsv", GTDBTK_CONTENT)
            result = self.mod.parse_gtdbtk(gtdbtk)
        self.assertFalse(self.mod.ENTEROBACTERALES not in result["bin_001"])
        self.assertTrue(self.mod.ENTEROBACTERALES not in result["bin_002"])

    def test_is_unexpected_taxon_unclassified(self):
        self.assertFalse(self.mod.is_unexpected_taxon("unclassified"))

    def test_is_unexpected_taxon_classified_other_order(self):
        self.assertTrue(self.mod.is_unexpected_taxon(
            "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;"
            "f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus acidophilus"
        ))

    def test_is_unexpected_taxon_enterobacterales(self):
        self.assertFalse(self.mod.is_unexpected_taxon(
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"
        ))

    def test_end_to_end_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            checkm2 = self._write(tmp, "quality.tsv", CHECKM2_CONTENT)
            gtdbtk = self._write(tmp, "gtdbtk.tsv", GTDBTK_CONTENT)
            tblout_dir = os.path.join(tmp, "tblout")
            os.makedirs(tblout_dir)
            self._write(tblout_dir, "bin_001.tblout", TBLOUT_PKS_CONTENT)
            context_dir = os.path.join(tmp, "context")
            os.makedirs(context_dir)
            self._write(context_dir, "bin_001.context.tsv", CONTEXT_CONTENT)
            out = os.path.join(tmp, "summary.tsv")
            import sys
            sys.argv = ["build_mag_summary.py",
                        "--checkm2", checkm2, "--gtdbtk", gtdbtk,
                        "--tblout_dir", tblout_dir, "--context_dir", context_dir,
                        "--sample", "SAMPLE1", "--evalue", "1e-5", "--out", out]
            self.mod.main()
            rows = list(csv.DictReader(open(out), delimiter="\t"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sample"], "SAMPLE1")
        self.assertEqual(rows[0]["bin_id"], "bin_001")
        self.assertEqual(rows[0]["clb_genes_detected"], "2")
        self.assertIn("clbA", rows[0]["clb_genes"])
        self.assertIn("clbB", rows[0]["clb_genes"])
        self.assertEqual(rows[0]["has_integrase"], "True")
        self.assertEqual(rows[0]["unexpected_taxon_flag"], "False")


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from molecular_parser.analysis import analyze_fasta, analyze_record, write_csv


class AnalysisTests(unittest.TestCase):
    def test_auto_detects_dna_and_counts_bases(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            fasta_path = Path(temp_dir, "dna.fasta")
            fasta_path.write_text(">dna_sample\nATGCGC\n", encoding="utf-8")

            reports = analyze_fasta(fasta_path)

        report = reports[0]
        self.assertEqual(report.sequence_type, "dna")
        self.assertEqual(report.length, 6)
        self.assertEqual(report.g_count, 2)
        self.assertEqual(report.c_count, 2)
        self.assertAlmostEqual(report.gc_content or 0.0, 66.6667, places=3)

    def test_forced_rna_analysis_uses_u_counts(self) -> None:
        report = analyze_record(
            SeqRecord(Seq("AUGCUU"), id="rna_sample", description="rna_sample"),
            sequence_type="rna",
        )

        self.assertEqual(report.sequence_type, "rna")
        self.assertEqual(report.u_count, 3)
        self.assertEqual(report.t_count, 0)
        self.assertIsNotNone(report.molecular_weight)

    def test_protein_analysis_populates_protein_metrics(self) -> None:
        report = analyze_record(
            SeqRecord(
                Seq("VLSPADKTNVKAAW"),
                id="protein_sample",
                description="protein_sample",
            ),
            sequence_type="protein",
        )

        self.assertEqual(report.sequence_type, "protein")
        self.assertIsNone(report.gc_content)
        self.assertIsNotNone(report.isoelectric_point)
        self.assertIsNotNone(report.gravy)
        self.assertIsNotNone(report.molecular_weight)

    def test_write_csv_exports_expected_headers(self) -> None:
        report = analyze_record(
            SeqRecord(Seq("ATGC"), id="dna_sample", description="dna_sample"),
            sequence_type="dna",
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            csv_path = write_csv([report], Path(temp_dir, "report.csv"))
            with csv_path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["record_id"], "dna_sample")
        self.assertEqual(rows[0]["sequence_type"], "dna")


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from molecular_parser import cli


class CliTests(unittest.TestCase):
    def test_resolve_input_files_filters_supported_extensions(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            fasta_file = root / "sample.fasta"
            fasta_file.write_text(">sample\nATGC\n", encoding="utf-8")
            (root / "notes.txt").write_text("ignore", encoding="utf-8")

            result = cli._resolve_input_files(root)

        self.assertEqual(result, [fasta_file])

    def test_main_uses_default_input_and_output_folders(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_dir = root / "input"
            output_dir = root / "output"
            fasta_file = input_dir / "sample.fasta"
            input_dir.mkdir()
            fasta_file.write_text(">sample\nATGC\n", encoding="utf-8")

            report = SimpleNamespace(
                record_id="sample",
                sequence_type="dna",
                length=4,
            )

            previous_cwd = Path.cwd()
            os.chdir(root)
            try:
                with patch("molecular_parser.cli.analyze_fasta", return_value=[report]) as mock_analyze:
                    with patch("molecular_parser.cli.write_csv", side_effect=lambda reports, path: path) as mock_write:
                        exit_code = cli.main([])
            finally:
                os.chdir(previous_cwd)

        self.assertEqual(exit_code, 0)
        mock_analyze.assert_called_once_with(fasta_file, sequence_type="auto")
        mock_write.assert_called_once_with([report], output_dir / "sample.csv")

    def test_main_creates_default_folders_when_missing(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            previous_cwd = Path.cwd()
            os.chdir(root)
            try:
                exit_code = cli.main([])
                input_exists = (root / "input").exists()
                output_exists = (root / "output").exists()
            finally:
                os.chdir(previous_cwd)

        self.assertEqual(exit_code, 1)
        self.assertTrue(input_exists)
        self.assertTrue(output_exists)


if __name__ == "__main__":
    unittest.main()

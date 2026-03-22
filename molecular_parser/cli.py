from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .analysis import SequenceType, analyze_fasta, write_csv

DEFAULT_INPUT_DIR = Path("input")
DEFAULT_OUTPUT_DIR = Path("output")
FASTA_SUFFIXES = {".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn", ".fas"}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze FASTA DNA/RNA/protein sequences and optionally export CSV.",
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Path to a FASTA file or folder. Defaults to input/.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Optional CSV output file for a single FASTA input.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Folder for generated CSV reports. Defaults to output/.",
    )
    parser.add_argument(
        "--sequence-type",
        choices=("auto", "dna", "rna", "protein"),
        default="auto",
        help="Override auto-detection when the FASTA content is ambiguous.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    input_path = args.input

    if not input_path.exists():
        if input_path == DEFAULT_INPUT_DIR:
            _ensure_default_folders()
            print("Created input/ and output/. Put FASTA files in input/ and run again.")
            return 1
        print(f"Input path not found: {input_path}")
        return 1

    try:
        input_files = _resolve_input_files(input_path)
    except ValueError as exc:
        print(str(exc))
        return 1

    if not input_files:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        print(
            f"No FASTA files found in {input_path}. "
            f"Supported extensions: {', '.join(sorted(FASTA_SUFFIXES))}"
        )
        return 1

    if args.output and len(input_files) != 1:
        print("--output can only be used with a single FASTA file.")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for fasta_path in input_files:
        try:
            reports = analyze_fasta(fasta_path, sequence_type=args.sequence_type)
        except ValueError as exc:
            print(f"{fasta_path}: {exc}")
            return 1

        _print_summary(fasta_path, reports, args.sequence_type)
        destination = _resolve_output_path(fasta_path, args.output, args.output_dir)
        destination = write_csv(reports, destination)
        print(f"CSV report written to {destination}")

    return 0


def _ensure_default_folders() -> None:
    DEFAULT_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _resolve_input_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        if input_path.suffix.lower() not in FASTA_SUFFIXES:
            raise ValueError(
                f"Unsupported input file extension: {input_path.suffix or '<none>'}"
            )
        return [input_path]

    if not input_path.is_dir():
        raise ValueError(f"Input path is neither a file nor a directory: {input_path}")

    return sorted(
        path
        for path in input_path.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )


def _resolve_output_path(
    input_path: Path,
    explicit_output: Path | None,
    output_dir: Path,
) -> Path:
    if explicit_output is not None:
        return explicit_output
    return output_dir / f"{input_path.stem}.csv"


def _print_summary(input_path: Path, reports: list, sequence_type: SequenceType) -> None:
    print(f"Analyzed {len(reports)} sequence(s) from {input_path}")
    print(f"Sequence type mode: {sequence_type}")
    for report in reports:
        print(
            f"- {report.record_id}: {report.sequence_type} | length={report.length}"
        )

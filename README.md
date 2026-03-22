# molecular-parser

`molecular-parser` is a fresh rebuild of my old Biopython-based project. It reads FASTA files, detects DNA/RNA/protein sequences, calculates useful sequence metrics, and exports the results as CSV.

## Features

- Parse multi-record FASTA files with Biopython `SeqIO`
- Auto-detect DNA, RNA, or protein sequences, with a CLI override when needed
- Report sequence length, molecular weight, GC content, nucleotide counts, and protein properties
- Export analysis results to CSV for spreadsheet workflows

## Install

For the simplest local setup:

```bash
python -m pip install biopython
```

If you also want the `molecular-parser` console command:

```bash
python -m pip install -e .
```

## Quick Start

```bash
python -m molecular_parser
```

The project now includes two fixed folders:

- `input/` for FASTA files
- `output/` for generated CSV reports

Drop your FASTA files into `input/` and run the command above from the project root. A sample file is already included at [input/sample_sequences.fasta](input/sample_sequences.fasta).

```fasta
>dna_example Human beta-globin DNA fragment
ATGGTGCACCTGACTCCTGAGGAGAAGTCT
>rna_example Example RNA fragment
AUGGUGCACCUGACUCCUGAGGAGAAGUCU
>protein_example Example peptide
VLSPADKTNVKAAWGKVGAHAGEYGAEALER
```

Running `python -m molecular_parser` will read all supported FASTA files from `input/` and create matching CSV files in `output/`.

Supported FASTA extensions:

- `.fasta`
- `.fa`
- `.fna`
- `.ffn`
- `.faa`
- `.frn`
- `.fas`

## Custom Paths

You can still point the parser to a specific FASTA file or a different folder:

```bash
python -m molecular_parser input/sample_sequences.fasta
```

Force a specific sequence type when auto-detection is too ambiguous:

```bash
python -m molecular_parser input/sample_sequences.fasta --sequence-type protein
```

Use a different output folder:

```bash
python -m molecular_parser path/to/fasta-folder --output-dir reports
```

Use an explicit CSV filename for a single FASTA file:

```bash
python -m molecular_parser input/sample_sequences.fasta --output reports/custom.csv
```

Run the test suite with:

```bash
python -m unittest discover -s tests
```

## Notes

- Auto-detection is heuristic. Sequences composed only of characters shared by nucleotide ambiguity codes and amino-acid symbols may need `--sequence-type`.
- Biopython is required for runtime and tests.

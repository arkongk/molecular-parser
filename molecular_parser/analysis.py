from __future__ import annotations

import csv
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Literal

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

DetectedSequenceType = Literal["dna", "rna", "protein"]
SequenceType = Literal["auto", "dna", "rna", "protein"]

NUCLEOTIDE_LETTERS = set("ACGTUNRYKMSWBDHV")
PROTEIN_ANALYSIS_LETTERS = set("ACDEFGHIKLMNPQRSTVWY")
CSV_FIELDNAMES = [
    "record_id",
    "description",
    "sequence_type",
    "length",
    "gc_content",
    "molecular_weight",
    "a_count",
    "c_count",
    "g_count",
    "t_count",
    "u_count",
    "n_count",
    "isoelectric_point",
    "aromaticity",
    "instability_index",
    "gravy",
    "helix_fraction",
    "turn_fraction",
    "sheet_fraction",
    "notes",
]


@dataclass(slots=True)
class SequenceReport:
    record_id: str
    description: str
    sequence_type: DetectedSequenceType
    length: int
    gc_content: float | None = None
    molecular_weight: float | None = None
    a_count: int | None = None
    c_count: int | None = None
    g_count: int | None = None
    t_count: int | None = None
    u_count: int | None = None
    n_count: int | None = None
    isoelectric_point: float | None = None
    aromaticity: float | None = None
    instability_index: float | None = None
    gravy: float | None = None
    helix_fraction: float | None = None
    turn_fraction: float | None = None
    sheet_fraction: float | None = None
    notes: str = ""

    def to_row(self) -> dict[str, object]:
        row = asdict(self)
        for key, value in row.items():
            if isinstance(value, float):
                row[key] = round(value, 4)
        return row


def analyze_fasta(path: str | Path, sequence_type: SequenceType = "auto") -> list[SequenceReport]:
    fasta_path = Path(path)
    with fasta_path.open("r", encoding="utf-8") as handle:
        return [
            analyze_record(record, sequence_type=sequence_type)
            for record in SeqIO.parse(handle, "fasta")
        ]


def analyze_record(record: SeqRecord, sequence_type: SequenceType = "auto") -> SequenceReport:
    sequence, cleanup_notes = _sanitize_sequence(str(record.seq))
    if not sequence:
        raise ValueError(f"Record '{record.id}' does not contain a sequence.")

    resolved_type = detect_sequence_type(sequence) if sequence_type == "auto" else sequence_type
    if resolved_type in {"dna", "rna"}:
        report = _analyze_nucleic_acid(record, sequence, resolved_type)
    else:
        report = _analyze_protein(record, sequence)

    if cleanup_notes:
        notes = cleanup_notes[:]
        if report.notes:
            notes.append(report.notes)
        report.notes = "; ".join(notes)
    return report


def detect_sequence_type(sequence: str) -> DetectedSequenceType:
    letters = {character for character in sequence.upper() if character.isalpha()}
    if not letters:
        raise ValueError("Cannot detect sequence type from an empty sequence.")
    if letters <= NUCLEOTIDE_LETTERS:
        return "rna" if "U" in letters and "T" not in letters else "dna"
    return "protein"


def write_csv(reports: Iterable[SequenceReport], output_path: str | Path) -> Path:
    rows = list(reports)
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        for report in rows:
            writer.writerow(report.to_row())

    return destination


def _analyze_nucleic_acid(
    record: SeqRecord,
    sequence: str,
    sequence_type: Literal["dna", "rna"],
) -> SequenceReport:
    counts = Counter(sequence)
    notes: list[str] = []

    try:
        weight = molecular_weight(sequence, seq_type=sequence_type)
    except ValueError as exc:
        weight = None
        notes.append(str(exc))

    return SequenceReport(
        record_id=record.id,
        description=record.description,
        sequence_type=sequence_type,
        length=len(sequence),
        gc_content=gc_fraction(sequence) * 100,
        molecular_weight=weight,
        a_count=counts.get("A", 0),
        c_count=counts.get("C", 0),
        g_count=counts.get("G", 0),
        t_count=counts.get("T", 0),
        u_count=counts.get("U", 0),
        n_count=counts.get("N", 0),
        notes="; ".join(notes),
    )


def _analyze_protein(record: SeqRecord, sequence: str) -> SequenceReport:
    notes: list[str] = []
    residue_sequence = sequence.replace("*", "")
    if residue_sequence != sequence:
        notes.append("Stop codons were excluded from protein metrics.")

    if not residue_sequence:
        return SequenceReport(
            record_id=record.id,
            description=record.description,
            sequence_type="protein",
            length=0,
            notes="; ".join(notes + ["Protein sequence is empty after cleanup."]),
        )

    unsupported = sorted({character for character in residue_sequence if character not in PROTEIN_ANALYSIS_LETTERS})
    protein_metrics: dict[str, float | None] = {
        "isoelectric_point": None,
        "aromaticity": None,
        "instability_index": None,
        "gravy": None,
        "helix_fraction": None,
        "turn_fraction": None,
        "sheet_fraction": None,
    }

    if unsupported:
        notes.append(
            "Protein metrics skipped for unsupported residues: "
            + "".join(unsupported)
        )
    else:
        protein = ProteinAnalysis(residue_sequence)
        helix, turn, sheet = protein.secondary_structure_fraction()
        protein_metrics = {
            "isoelectric_point": protein.isoelectric_point(),
            "aromaticity": protein.aromaticity(),
            "instability_index": protein.instability_index(),
            "gravy": protein.gravy(),
            "helix_fraction": helix,
            "turn_fraction": turn,
            "sheet_fraction": sheet,
        }

    try:
        weight = molecular_weight(residue_sequence, seq_type="protein")
    except ValueError as exc:
        weight = None
        notes.append(str(exc))

    return SequenceReport(
        record_id=record.id,
        description=record.description,
        sequence_type="protein",
        length=len(residue_sequence),
        molecular_weight=weight,
        notes="; ".join(notes),
        **protein_metrics,
    )


def _sanitize_sequence(raw_sequence: str) -> tuple[str, list[str]]:
    sequence = raw_sequence.upper()
    notes: list[str] = []

    if any(character.isspace() for character in sequence):
        sequence = "".join(character for character in sequence if not character.isspace())

    gap_characters = {"-", "."}
    if any(character in gap_characters for character in sequence):
        sequence = "".join(character for character in sequence if character not in gap_characters)
        notes.append("Gap characters were removed before analysis.")

    return sequence, notes

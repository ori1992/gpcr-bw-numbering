"""GPCR BW Annotator: annotate a class-A GPCR sequence and export a table.

The application is deliberately self-contained.  It uses only the Python
standard library, performs a BLOSUM62/Gotoh pairwise alignment to bovine
rhodopsin, transfers a class-A topology, validates the canonical BW anchors,
and writes portable HTML and tab-separated reports without network access.
"""

from __future__ import annotations

import argparse
import re
import sys
import threading
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, Sequence
from xml.sax.saxutils import escape


APP_NAME = "GPCR BW Annotator"
APP_VERSION = "1.3.0"

BOVINE_HEADER = "sp|P02699|OPSD_BOVIN Rhodopsin OS=Bos taurus"
BOVINE_SEQUENCE = (
    "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLY"
    "VTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLG"
    "GEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIP"
    "EGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQES"
    "ATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAV"
    "YNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA"
)

# Curated GPCRdb opsd_bovin boundaries (official residue service, 2026-08-18).
# They are transferred through the pairwise alignment; they are not copied as
# target residue numbers.
BOVINE_SEGMENTS: tuple[tuple[str, int, int], ...] = (
    ("N-terminus", 1, 32),
    ("TM1", 33, 65),
    ("ICL1", 66, 69),
    ("TM2", 70, 101),
    ("ECL1", 102, 105),
    ("TM3", 106, 141),
    ("ICL2", 142, 148),
    ("TM4", 149, 174),
    ("ECL2", 175, 198),
    ("TM5", 199, 237),
    ("ICL3", 238, 239),
    ("TM6", 240, 278),
    ("ECL3", 279, 283),
    ("TM7", 284, 309),
    ("Helix 8", 310, 322),
    ("C-terminus", 323, 348),
)

LOOP_CODES = {
    "ICL1": "12",
    "ECL1": "23",
    "ICL2": "34",
    "ECL2": "45",
    "ICL3": "56",
    "ECL3": "67",
}

# ICL1/ECL1/ICL2/ECL2 anchors preserve curated bovine GPCRdb positions.
# GPCRdb does not number every flexible ICL3/ECL3 residue; their anchors define
# this application's explicitly alignment-relative continuation only.
LOOP_TEMPLATE_ANCHORS = {
    "ICL1": 68,
    "ECL1": 103,
    "ICL2": 142,
    "ECL2": 187,
    "ICL3": 239,
    "ECL3": 281,
}

ANCHORS: tuple[tuple[str, int, str, str], ...] = (
    ("1.50", 55, "N", "TM1"),
    ("2.50", 83, "D", "TM2"),
    ("3.50", 135, "R", "TM3"),
    ("4.50", 161, "W", "TM4"),
    ("45.50", 187, "C", "ECL2"),
    ("5.50", 215, "P", "TM5"),
    ("6.50", 267, "P", "TM6"),
    ("7.50", 303, "P", "TM7"),
)

HYDROPHOBIC = set("ACFILMVWY")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")

BLOSUM62_TEXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1
"""


def _parse_blosum62() -> dict[tuple[str, str], int]:
    lines = [line.split() for line in BLOSUM62_TEXT.strip().splitlines()]
    alphabet = lines[0]
    scores: dict[tuple[str, str], int] = {}
    for row in lines[1:]:
        aa = row[0]
        for other, value in zip(alphabet, row[1:]):
            scores[(aa, other)] = int(value)
    return scores


BLOSUM62 = _parse_blosum62()


@dataclass(frozen=True)
class SequenceInput:
    name: str
    sequence: str


@dataclass(frozen=True)
class AnchorResult:
    label: str
    position: int
    residue: str
    expected: str
    segment: str
    matched: bool
    source: str = "automatic"


@dataclass(frozen=True)
class ResidueAnnotation:
    position: int
    residue: str
    segment: str
    generic_number: str
    annotation: str


@dataclass(frozen=True)
class TopologyBoundaries:
    automatic_tm5_end: int
    automatic_tm6_start: int
    tm5_end: int
    tm6_start: int
    tm5_end_source: str = "automatic"
    tm6_start_source: str = "automatic"

    @property
    def manually_adjusted(self) -> bool:
        return self.tm5_end_source == "manual" or self.tm6_start_source == "manual"


@dataclass(frozen=True)
class AnalysisResult:
    sequence_name: str
    sequence: str
    annotations: tuple[ResidueAnnotation, ...]
    anchors: tuple[AnchorResult, ...]
    confidence: str
    confidence_detail: str
    identity: float
    tm_hydrophobicity: float
    warnings: tuple[str, ...]
    aligned_template: str
    aligned_target: str
    topology_boundaries: TopologyBoundaries


def parse_sequence_text(text: str, fallback_name: str = "GPCR_sequence") -> SequenceInput:
    raw = text.strip()
    if not raw:
        raise ValueError("Paste a protein sequence or load a FASTA file first.")
    headers = [line.strip()[1:].strip() for line in raw.splitlines() if line.strip().startswith(">")]
    if len(headers) > 1:
        raise ValueError("This focused annotator accepts one protein sequence at a time.")
    name = headers[0] if headers else fallback_name
    sequence_lines = [line.strip() for line in raw.splitlines() if not line.strip().startswith(">")]
    sequence = re.sub(r"[\s\d\-\.]+", "", "".join(sequence_lines)).upper().rstrip("*")
    invalid = sorted(set(sequence) - VALID_AA)
    if invalid:
        raise ValueError("Unsupported sequence character(s): " + ", ".join(invalid))
    if len(sequence) < 180:
        raise ValueError("The sequence is too short for a complete class-A seven-helix GPCR core (minimum 180 aa).")
    if len(sequence) > 1200:
        raise ValueError("The sequence is longer than 1,200 aa; provide the GPCR protein chain only.")
    return SequenceInput(name=name or fallback_name, sequence=sequence)


def read_sequence_file(path: Path) -> SequenceInput:
    try:
        text = path.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError:
        text = path.read_text(encoding="latin-1")
    return parse_sequence_text(text, fallback_name=path.stem)


def pairwise_align(template: str, target: str) -> tuple[str, str, int]:
    """Global Gotoh alignment with BLOSUM62 and affine gap penalties."""

    n, m = len(template), len(target)
    neg = -10**9
    gap_open, gap_extend = 10, 1
    mat = [[neg] * (m + 1) for _ in range(n + 1)]
    gap_target = [[neg] * (m + 1) for _ in range(n + 1)]
    gap_template = [[neg] * (m + 1) for _ in range(n + 1)]
    trace_m = [bytearray(m + 1) for _ in range(n + 1)]
    trace_x = [bytearray(m + 1) for _ in range(n + 1)]
    trace_y = [bytearray(m + 1) for _ in range(n + 1)]
    mat[0][0] = 0
    for i in range(1, n + 1):
        gap_target[i][0] = -gap_open - (i - 1) * gap_extend
        trace_x[i][0] = 1 if i > 1 else 0
    for j in range(1, m + 1):
        gap_template[0][j] = -gap_open - (j - 1) * gap_extend
        trace_y[0][j] = 2 if j > 1 else 0

    for i in range(1, n + 1):
        aa = template[i - 1]
        mat_i, x_i, y_i = mat[i], gap_target[i], gap_template[i]
        mat_prev, x_prev, y_prev = mat[i - 1], gap_target[i - 1], gap_template[i - 1]
        tm_i, tx_i, ty_i = trace_m[i], trace_x[i], trace_y[i]
        for j in range(1, m + 1):
            previous = (mat_prev[j - 1], x_prev[j - 1], y_prev[j - 1])
            state = max(range(3), key=previous.__getitem__)
            score = BLOSUM62.get((aa, target[j - 1]), -1)
            mat_i[j] = previous[state] + score
            tm_i[j] = state

            candidates_x = (mat_prev[j] - gap_open, x_prev[j] - gap_extend, y_prev[j] - gap_open)
            state_x = max(range(3), key=candidates_x.__getitem__)
            x_i[j] = candidates_x[state_x]
            tx_i[j] = state_x

            candidates_y = (mat_i[j - 1] - gap_open, x_i[j - 1] - gap_open, y_i[j - 1] - gap_extend)
            state_y = max(range(3), key=candidates_y.__getitem__)
            y_i[j] = candidates_y[state_y]
            ty_i[j] = state_y

    final_scores = (mat[n][m], gap_target[n][m], gap_template[n][m])
    state = max(range(3), key=final_scores.__getitem__)
    score = final_scores[state]
    i, j = n, m
    aligned_template: list[str] = []
    aligned_target: list[str] = []
    while i > 0 or j > 0:
        if state == 0:
            if i == 0:
                state = 2
                continue
            if j == 0:
                state = 1
                continue
            previous_state = trace_m[i][j]
            aligned_template.append(template[i - 1])
            aligned_target.append(target[j - 1])
            i -= 1
            j -= 1
        elif state == 1:
            if i == 0:
                state = 2
                continue
            previous_state = trace_x[i][j]
            aligned_template.append(template[i - 1])
            aligned_target.append("-")
            i -= 1
        else:
            if j == 0:
                state = 1
                continue
            previous_state = trace_y[i][j]
            aligned_template.append("-")
            aligned_target.append(target[j - 1])
            j -= 1
        state = previous_state
    return "".join(reversed(aligned_template)), "".join(reversed(aligned_target)), score


def _residue_column(aligned_sequence: str, residue_number: int) -> int:
    count = 0
    for column, character in enumerate(aligned_sequence):
        if character != "-":
            count += 1
            if count == residue_number:
                return column
    raise ValueError(f"Residue {residue_number} is outside the aligned sequence.")


def _target_position_at_column(aligned_target: str, column: int) -> int:
    if aligned_target[column] != "-":
        return sum(character != "-" for character in aligned_target[: column + 1])
    for radius in range(1, 80):
        for candidate in (column - radius, column + radius):
            if 0 <= candidate < len(aligned_target) and aligned_target[candidate] != "-":
                return sum(character != "-" for character in aligned_target[: candidate + 1])
    raise ValueError("Could not transfer a reference position through the alignment.")


def _anchor_context_score(label: str, target: str, position: int) -> float:
    """Score conserved class-A/opsin context at a 1-based anchor candidate."""

    expected = {item_label: identity for item_label, _pos, identity, _segment in ANCHORS}[label]
    index = position - 1
    score = 7.0 if target[index] == expected else 0.0
    if label == "3.50":
        if index >= 1 and target[index - 1] in "DE":
            score += 3.0
        if index + 1 < len(target) and target[index + 1] in "YFWH":
            score += 3.0
    elif label == "6.50":
        if index >= 2 and target[index - 2] == "W":
            score += 4.0
        if index >= 3 and target[index - 3] in "CF":
            score += 2.0
    elif label == "7.50":
        if index >= 1 and target[index - 1] == "N":
            score += 4.0
        if index + 3 < len(target) and target[index + 3] in "YF":
            score += 4.0
    elif label == "45.50":
        neighborhood = target[max(0, index - 4): min(len(target), index + 5)]
        if neighborhood.count("C") >= 2:
            score += 1.0
    return score


def _detect_anchors(aligned_template: str, aligned_target: str, target: str) -> list[AnchorResult]:
    results_by_label: dict[str, AnchorResult] = {}
    # Detect the seven helix anchors first.  ECL2 is handled separately because
    # many receptors contain a local cysteine cluster around the conserved
    # disulfide-forming cysteine.
    for label, template_position, expected, segment in ANCHORS:
        if label == "45.50":
            continue
        template_column = _residue_column(aligned_template, template_position)
        projected = _target_position_at_column(aligned_target, template_column)
        candidates = range(max(1, projected - 12), min(len(target), projected + 12) + 1)
        target_position = max(
            candidates,
            key=lambda position: (
                _anchor_context_score(label, target, position) - 0.45 * abs(position - projected),
                -abs(position - projected),
            ),
        )
        residue = target[target_position - 1]
        results_by_label[label] = AnchorResult(label, target_position, residue, expected, segment, residue == expected)

    # In bovine rhodopsin C45.50 lies 26 residues after W4.50 within the
    # 54-residue 4.50-to-5.50 interval.  Interpolating that relative location
    # in the target is much more reliable than selecting the alignment-nearest
    # cysteine when ECL2 contains adjacent or clustered cysteines (e.g. ADRB2).
    p4 = results_by_label["4.50"].position
    p5 = results_by_label["5.50"].position
    predicted_ecl2 = round(p4 + (p5 - p4) * ((187 - 161) / (215 - 161)))
    ecl2_candidates = [
        position
        for position in range(max(p4 + 1, predicted_ecl2 - 18), min(p5, predicted_ecl2 + 18) + 1)
        if target[position - 1] == "C"
    ]
    if ecl2_candidates:
        ecl2_column = _residue_column(aligned_template, 187)
        aligned_ecl2 = _target_position_at_column(aligned_target, ecl2_column)
        ecl2_position = min(
            ecl2_candidates,
            key=lambda position: (
                abs(position - predicted_ecl2),
                abs(position - aligned_ecl2),
                position,
            ),
        )
    else:
        ecl2_column = _residue_column(aligned_template, 187)
        ecl2_position = _target_position_at_column(aligned_target, ecl2_column)
    ecl2_residue = target[ecl2_position - 1]
    results_by_label["45.50"] = AnchorResult(
        "45.50", ecl2_position, ecl2_residue, "C", "ECL2", ecl2_residue == "C"
    )
    return [results_by_label[label] for label, _position, _expected, _segment in ANCHORS]


def parse_manual_anchors(text: str) -> dict[str, int]:
    if not text.strip():
        return {}
    allowed = {label for label, *_rest in ANCHORS}
    mapping: dict[str, int] = {}
    for piece in re.split(r"[,;\n]+", text):
        if not piece.strip():
            continue
        match = re.fullmatch(r"\s*(\d{1,2}\.\d{1,3})\s*(?:=|:)\s*(\d+)\s*", piece)
        if not match:
            raise ValueError(f"Could not read manual anchor '{piece.strip()}'. Use 3.50=135.")
        label = match.group(1)
        left, right = label.split(".", 1)
        label = f"{left}.{int(right):02d}"
        if label not in allowed:
            raise ValueError(f"{label} is not one of the supported .50/ECL2 anchors.")
        mapping[label] = int(match.group(2))
    return mapping


def _apply_manual_anchors(anchors: Sequence[AnchorResult], manual: dict[str, int], target: str) -> list[AnchorResult]:
    if not manual:
        return list(anchors)
    updated: list[AnchorResult] = []
    for item in anchors:
        if item.label not in manual:
            updated.append(item)
            continue
        position = manual[item.label]
        if not 1 <= position <= len(target):
            raise ValueError(f"Manual anchor {item.label}={position} is outside the sequence.")
        residue = target[position - 1]
        updated.append(
            AnchorResult(item.label, position, residue, item.expected, item.segment, residue == item.expected, "manual")
        )
    return updated


def _segment_for_template_residue(position: int) -> str:
    for name, start, end in BOVINE_SEGMENTS:
        if start <= position <= end:
            return name
    return "Unassigned"


def transfer_segments(aligned_template: str, aligned_target: str, target_length: int) -> list[str]:
    column_segments: list[str | None] = []
    template_position = 0
    for character in aligned_template:
        if character == "-":
            column_segments.append(None)
        else:
            template_position += 1
            column_segments.append(_segment_for_template_residue(template_position))

    left: list[str | None] = []
    current: str | None = None
    for item in column_segments:
        if item is not None:
            current = item
        left.append(current)
    right: list[str | None] = [None] * len(column_segments)
    current = None
    for index in range(len(column_segments) - 1, -1, -1):
        if column_segments[index] is not None:
            current = column_segments[index]
        right[index] = current

    def choose_gap_segment(before: str | None, after: str | None) -> str:
        if before is None:
            return after or "N-terminus"
        if after is None:
            return before
        if before == after:
            return before
        before_tm = before.startswith("TM")
        after_tm = after.startswith("TM")
        if before_tm != after_tm:
            return after if not after_tm else before
        return before

    result: list[str] = []
    for column, target_character in enumerate(aligned_target):
        if target_character == "-":
            continue
        segment = column_segments[column] or choose_gap_segment(left[column], right[column])
        result.append(segment)
    if len(result) != target_length:
        raise RuntimeError("Internal topology transfer length mismatch.")
    return result


def _insertion_suffix(index: int) -> str:
    value = index + 1
    pieces: list[str] = []
    while value:
        value, remainder = divmod(value - 1, 26)
        pieces.append(chr(ord("a") + remainder))
    return "".join(reversed(pieces))


def _template_loop_label(segment: str, template_position: int) -> str:
    index = 50 + template_position - LOOP_TEMPLATE_ANCHORS[segment]
    return f"{LOOP_CODES[segment]}.{index:02d}"


def assign_loop_numbers(
    aligned_template: str,
    aligned_target: str,
    effective_segments: Sequence[str],
    anchor_map: dict[str, int],
) -> dict[int, str]:
    """Number every target ECL/ICL residue without shifting homologous labels.

    Template-aligned residues retain the bovine-reference loop label. Insertions
    use alphabetic suffixes after the preceding aligned position. Residues moved
    into ICL3 by manual TM5/TM6 boundaries are numbered relative to the projected
    56.50 anchor, so the manual interval remains completely annotated.
    """

    target_loop_anchors: dict[str, int] = {}
    for segment, template_anchor in LOOP_TEMPLATE_ANCHORS.items():
        if segment == "ECL2":
            target_loop_anchors[segment] = anchor_map["45.50"]
        else:
            column = _residue_column(aligned_template, template_anchor)
            target_loop_anchors[segment] = _target_position_at_column(aligned_target, column)

    template_position = 0
    target_position = 0
    template_at_column: list[int | None] = []
    target_at_column: list[int | None] = []
    for template_aa, target_aa in zip(aligned_template, aligned_target):
        if template_aa != "-":
            template_position += 1
            template_at_column.append(template_position)
        else:
            template_at_column.append(None)
        if target_aa != "-":
            target_position += 1
            target_at_column.append(target_position)
        else:
            target_at_column.append(None)

    labels: dict[int, str] = {}
    insertion_counts: dict[str, int] = {}
    used: set[str] = set()
    for column, position in enumerate(target_at_column):
        if position is None:
            continue
        segment = effective_segments[position - 1]
        if segment not in LOOP_CODES:
            continue
        template_residue = template_at_column[column]
        if template_residue is not None and _segment_for_template_residue(template_residue) == segment:
            label = _template_loop_label(segment, template_residue)
        elif template_residue is None:
            previous = next(
                (
                    template_at_column[index]
                    for index in range(column - 1, -1, -1)
                    if template_at_column[index] is not None
                    and _segment_for_template_residue(int(template_at_column[index])) == segment
                ),
                None,
            )
            if previous is None:
                following = next(
                    (
                        template_at_column[index]
                        for index in range(column + 1, len(template_at_column))
                        if template_at_column[index] is not None
                        and _segment_for_template_residue(int(template_at_column[index])) == segment
                    ),
                    LOOP_TEMPLATE_ANCHORS[segment],
                )
                base_index = 49 + int(following) - LOOP_TEMPLATE_ANCHORS[segment]
                base_label = f"{LOOP_CODES[segment]}.{base_index:02d}"
            else:
                base_label = _template_loop_label(segment, int(previous))
            insertion_counts[base_label] = insertion_counts.get(base_label, 0) + 1
            label = f"{base_label}{_insertion_suffix(insertion_counts[base_label] - 1)}"
        else:
            index = 50 + position - target_loop_anchors[segment]
            label = f"{LOOP_CODES[segment]}.{index:02d}"
        if label in used:
            duplicate = 0
            base_label = label
            while label in used:
                label = f"{base_label}{_insertion_suffix(duplicate)}"
                duplicate += 1
        used.add(label)
        labels[position] = label

    # ECL2 is anchored biologically by the disulfide-forming C45.50.  In a
    # cysteine-rich loop (for example ADRB2), a global alignment can place the
    # bovine C45.50 column on a different nearby cysteine.  Shift the complete
    # alignment-relative ECL2 label set so the independently detected anchor is
    # exactly 45.50, while retaining insertion suffixes and relative spacing.
    ecl2_anchor_position = anchor_map["45.50"]
    current_anchor_label = labels.get(ecl2_anchor_position, "")
    anchor_match = re.fullmatch(r"45\.(\d+)([a-z]*)", current_anchor_label)
    if anchor_match:
        delta = 50 - int(anchor_match.group(1))
        if delta:
            for position, label in list(labels.items()):
                if effective_segments[position - 1] != "ECL2":
                    continue
                match = re.fullmatch(r"45\.(\d+)([a-z]*)", label)
                if match:
                    labels[position] = f"45.{int(match.group(1)) + delta:02d}{match.group(2)}"
    return labels


def parse_manual_boundary(value: str | int | None, label: str) -> int | None:
    """Parse an optional one-based sequence position from the desktop form or API."""

    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    text = str(value).strip()
    if not re.fullmatch(r"\d+", text):
        raise ValueError(f"{label} must be a whole-number sequence position, for example 248.")
    position = int(text)
    if position < 1:
        raise ValueError(f"{label} must be at least 1.")
    return position


def apply_tm5_tm6_boundaries(
    segments: Sequence[str],
    anchor_map: dict[str, int],
    manual_tm5_end: str | int | None = None,
    manual_tm6_start: str | int | None = None,
) -> tuple[list[str], TopologyBoundaries]:
    """Apply optional structure-informed cytoplasmic TM5/TM6 boundaries.

    The .50 anchors remain unchanged.  Only the interval from 5.50 through
    6.50 is repartitioned into TM5, ICL3, and TM6.
    """

    tm5_positions = [position for position, segment in enumerate(segments, start=1) if segment == "TM5"]
    tm6_positions = [position for position, segment in enumerate(segments, start=1) if segment == "TM6"]
    if not tm5_positions or not tm6_positions:
        raise ValueError("Could not infer automatic TM5/TM6 boundaries from the alignment.")

    automatic_tm5_end = max(tm5_positions)
    automatic_tm6_start = min(tm6_positions)
    parsed_tm5_end = parse_manual_boundary(manual_tm5_end, "TM5 end")
    parsed_tm6_start = parse_manual_boundary(manual_tm6_start, "TM6 start")
    tm5_end = parsed_tm5_end if parsed_tm5_end is not None else automatic_tm5_end
    tm6_start = parsed_tm6_start if parsed_tm6_start is not None else automatic_tm6_start
    anchor5 = anchor_map["5.50"]
    anchor6 = anchor_map["6.50"]

    if tm5_end < anchor5:
        raise ValueError(f"TM5 end ({tm5_end}) cannot be before the inferred 5.50 anchor ({anchor5}).")
    if tm6_start > anchor6:
        raise ValueError(f"TM6 start ({tm6_start}) cannot be after the inferred 6.50 anchor ({anchor6}).")
    if tm5_end >= tm6_start:
        raise ValueError(
            f"TM5 end ({tm5_end}) must be before TM6 start ({tm6_start}); "
            "the residues between them are assigned to ICL3."
        )

    updated = list(segments)
    for position in range(anchor5, anchor6 + 1):
        if position <= tm5_end:
            updated[position - 1] = "TM5"
        elif position < tm6_start:
            updated[position - 1] = "ICL3"
        else:
            updated[position - 1] = "TM6"

    boundaries = TopologyBoundaries(
        automatic_tm5_end=automatic_tm5_end,
        automatic_tm6_start=automatic_tm6_start,
        tm5_end=tm5_end,
        tm6_start=tm6_start,
        tm5_end_source="manual" if parsed_tm5_end is not None else "automatic",
        tm6_start_source="manual" if parsed_tm6_start is not None else "automatic",
    )
    return updated, boundaries


def _alignment_identity(aligned_template: str, aligned_target: str) -> float:
    paired = [(a, b) for a, b in zip(aligned_template, aligned_target) if a != "-" and b != "-"]
    return sum(a == b for a, b in paired) / len(paired) if paired else 0.0


def _core_alignment_identity(aligned_template: str, aligned_target: str) -> float:
    template_position = 0
    paired: list[tuple[str, str]] = []
    for template_aa, target_aa in zip(aligned_template, aligned_target):
        if template_aa != "-":
            template_position += 1
        if 33 <= template_position <= 322 and template_aa != "-" and target_aa != "-":
            paired.append((template_aa, target_aa))
    return sum(left == right for left, right in paired) / len(paired) if paired else 0.0


def analyze_sequence(
    sequence_input: SequenceInput,
    manual_anchor_text: str = "",
    status: Callable[[str], None] | None = None,
    manual_tm5_end: str | int | None = None,
    manual_tm6_start: str | int | None = None,
) -> AnalysisResult:
    status = status or (lambda _message: None)
    status("Aligning the sequence to the class-A bovine rhodopsin reference...")
    aligned_template, aligned_target, _score = pairwise_align(BOVINE_SEQUENCE, sequence_input.sequence)
    status("Transferring topology and validating conserved anchors...")
    anchors = _detect_anchors(aligned_template, aligned_target, sequence_input.sequence)
    manual = parse_manual_anchors(manual_anchor_text)
    anchors = _apply_manual_anchors(anchors, manual, sequence_input.sequence)
    automatic_segments = transfer_segments(aligned_template, aligned_target, len(sequence_input.sequence))

    anchor_map = {item.label: item.position for item in anchors}
    segments, topology_boundaries = apply_tm5_tm6_boundaries(
        automatic_segments,
        anchor_map,
        manual_tm5_end,
        manual_tm6_start,
    )
    helix_anchor_map = {f"TM{helix}": anchor_map[f"{helix}.50"] for helix in range(1, 8)}
    loop_numbers = assign_loop_numbers(aligned_template, aligned_target, segments, anchor_map)
    annotations: list[ResidueAnnotation] = []
    for position, (residue, segment) in enumerate(zip(sequence_input.sequence, segments), start=1):
        generic = ""
        detail = segment
        if segment in helix_anchor_map:
            helix = segment[-1]
            index = 50 + position - helix_anchor_map[segment]
            generic = f"{helix}.{index:02d}"
            if topology_boundaries.manually_adjusted and segment != automatic_segments[position - 1]:
                detail = "Ballesteros-Weinstein (manual helix boundary)"
            else:
                detail = "Ballesteros-Weinstein"
        elif re.fullmatch(r"[EI]CL[1-3]", segment):
            generic = loop_numbers[position]
            detail = f"Alignment-relative {segment} numbering (segment {LOOP_CODES[segment]})"
        elif segment == "Helix 8":
            generic = f"8.{47 + position - min(
                index for index, name in enumerate(segments, start=1) if name == 'Helix 8'
            ):02d}"
            detail = "Helix 8 alignment-relative numbering"
        annotations.append(ResidueAnnotation(position, residue, segment, generic, detail))

    identity = _core_alignment_identity(aligned_template, aligned_target)
    # Assess the automatically transferred membrane core so that a curated,
    # less-hydrophobic cytoplasmic helix extension does not lower confidence.
    tm_residues = [
        residue
        for residue, segment in zip(sequence_input.sequence, automatic_segments)
        if segment.startswith("TM")
    ]
    tm_hydrophobicity = sum(aa in HYDROPHOBIC for aa in tm_residues) / len(tm_residues) if tm_residues else 0.0
    matched = sum(item.matched for item in anchors)
    motif_supported = sum(
        _anchor_context_score(label, sequence_input.sequence, anchor_map[label]) >= 10
        for label in ("3.50", "6.50", "7.50")
    )
    positions = [item.position for item in anchors]
    ordered = all(right > left for left, right in zip(positions, positions[1:]))
    warnings: list[str] = []
    if not ordered:
        warnings.append("The inferred anchors are not in the expected TM1-to-TM7 sequence order.")
    if matched < 8:
        missing = ", ".join(item.label for item in anchors if not item.matched)
        warnings.append(f"Conserved anchor identity was not recovered at: {missing}.")
    if identity < 0.25:
        warnings.append("7TM-core identity to bovine rhodopsin is low; segment boundaries should be reviewed.")
    if motif_supported < 2:
        warnings.append("Fewer than two diagnostic DRY-like, CWxP, and NPxxY motif contexts were recovered.")
    if tm_hydrophobicity < 0.55:
        warnings.append("Transferred TM segments have unexpectedly low hydrophobic content.")
    if not 240 <= len(sequence_input.sequence) <= 700:
        warnings.append("Sequence length is outside the usual compact class-A GPCR range (240-700 aa).")
    if manual:
        warnings.append("One or more anchors were manually overridden; verify those positions independently.")
    if topology_boundaries.manually_adjusted:
        warnings.append(
            "Manual TM5/TM6 boundaries were applied; verify them against a structure, model, or curated alignment."
        )

    if matched >= 7 and motif_supported >= 2 and ordered and identity >= 0.28 and tm_hydrophobicity >= 0.58:
        confidence = "High"
    elif matched >= 6 and motif_supported >= 1 and ordered and identity >= 0.18 and tm_hydrophobicity >= 0.52:
        confidence = "Moderate"
    else:
        confidence = "Low"
    confidence_detail = (
        f"{matched}/8 conserved anchor identities; {motif_supported}/3 diagnostic motif contexts; "
        f"{identity:.1%} 7TM-core identity to P02699; "
        f"{tm_hydrophobicity:.1%} hydrophobic residues across transferred TM segments."
    )
    status("Annotation complete.")
    return AnalysisResult(
        sequence_name=sequence_input.name,
        sequence=sequence_input.sequence,
        annotations=tuple(annotations),
        anchors=tuple(anchors),
        confidence=confidence,
        confidence_detail=confidence_detail,
        identity=identity,
        tm_hydrophobicity=tm_hydrophobicity,
        warnings=tuple(warnings),
        aligned_template=aligned_template,
        aligned_target=aligned_target,
        topology_boundaries=topology_boundaries,
    )


---
name: gpcr-bw-numbering
description: Assign class-A GPCR BW numbers and complete ECL1–3, ICL1–3, and helix-8 numbering from pasted sequences, FASTA, PDB/mmCIF, PDB accessions, or protein links. Return residue tables using GPCRBWAnnotator v1.3 conventions, including ECL2 45.xx labels and insertion suffixes, with curated GPCRdb values kept separate. For other GPCR classes, identify the appropriate scheme instead of mislabeling it as BW.
---

# GPCR BW numbering

Produce a complete residue table including every numbering type in GPCRBWAnnotator v1.3, with a brief account of the receptor, chain, numbering source, confidence, and unresolved assignments. Use the bundled deterministic helper; its local scientific engine is extracted unchanged from that software.

## Numbering included by default

Always include the software's full `software_number` column in the delivered table and chat preview. Do not restrict the result to the seven TM anchors or to positions with curated GPCRdb numbers.

| Region | Software labels |
|---|---|
| TM1–TM7 | Sequence-based BW, e.g. 3.50, 6.48, 7.53 |
| ECL2 | Complete 45.xx series around C45.50, including loop insertions |
| ECL1 / ECL3 | 23.xx / 67.xx |
| ICL1 / ICL2 / ICL3 | 12.xx / 34.xx / 56.xx |
| Helix 8 | 8.xx, starting at 8.47 in the software's transferred H8 segment |
| Insertions in loops | Alphabetic suffixes such as 45.43a and 45.43b, preserving downstream labels |

Use `software_segment` to interpret each software label. These full-loop and H8 values are explicitly alignment-relative; keep their provenance visible. Termini and unsupported regions remain in the table with blank numbers. Supported TM5/TM6 boundary adjustments still number the entire resulting ICL3 without moving 5.50/6.50.

## Input and receptor selection

Accept a pasted amino-acid sequence, FASTA, PDB or mmCIF (including gzip), PDB accession, or protein link. Resolve RCSB, UniProt, and GPCRdb entry links to their sequence/structure data. For other websites, read the page and obtain its explicit sequence, FASTA download, PDB accession, or linked structure; never treat stripped page prose as protein sequence. Keep source URLs and identifiers.

Identify the GPCR class before treating an assignment as BW. BW is the native sequence-based scheme for class A; classes B, C, and F have different native schemes. If the class is unknown, investigate identity from the supplied source and curated databases. Motifs alone do not establish class or receptor identity.

For multiple records or chains, run `--list`, identify the receptor using structure metadata, and pass `--chain` or `--record`. Do not concatenate chains or choose the longest chain automatically. If two receptor chains remain plausible, report the candidates and ask for the intended chain while completing independent preparation.

## Assignment

Run `scripts/annotate_gpcr.py` with Python 3.10 or later. It needs only the standard library. All paths below are relative to this skill's folder; replace them with absolute paths when running from another working directory. Use an available Python 3.10+ interpreter. In Codex desktop, `load_workspace_dependencies` can locate its bundled runtime; other environments can use their installed Python.

```text
python scripts/annotate_gpcr.py --input receptor.fasta --out outputs/receptor_bw
python scripts/annotate_gpcr.py --sequence "AMINO_ACID_SEQUENCE" --out outputs/receptor_bw
python scripts/annotate_gpcr.py --input receptor.pdb --list
python scripts/annotate_gpcr.py --input receptor.pdb --chain R --out outputs/receptor_bw
python scripts/annotate_gpcr.py --pdb-id 6OYA --chain R --out outputs/6OYA_R_bw
python scripts/annotate_gpcr.py --url https://www.uniprot.org/uniprotkb/P07550/entry --gpcrdb-entry adrb2_human --out outputs/adrb2_bw
```

Prefer a curated GPCRdb reference when the receptor and species can be verified. Supply its actual entry name using `--gpcrdb-entry`; do not guess a similarly named receptor. The helper checks class A, endpoint consistency, alignment identity, and TM coverage, and keeps curated BW and GPCRdb numbers separate from the full `software_number` series. The latter remains included even when the curated reference leaves flexible loops blank. Substitutions are noted; insertions without a curated equivalent remain unassigned in the curated columns. Software loop labels do not establish structural equivalence for those insertions. Inspect `software_number_status` and `summary.software_numbering`, including any disagreement with curated anchors or boundaries. For remote access, offline caches, and complex constructs, read [references/input-and-output.md](references/input-and-output.md).

When a suitable curated reference is unavailable, use the bundled v1.3 inference. It aligns to bovine rhodopsin with BLOSUM62/Gotoh and validates anchors and motifs. Inspect confidence, anchor identities and order, motif support, and transferred topology before relying on the table. For anchor values, exact conventions, and limitations, read [references/numbering-rules.md](references/numbering-rules.md). A high heuristic score is supporting evidence, not independent validation of class or every residue's equivalence.

Preserve author residue numbers, insertion codes, chain IDs, and sequence positions as separate fields. PDB SEQRES and mmCIF polymer sequences retain unresolved positions. Observed-only structures need a matching full sequence through `--reference-fasta`; compacting ATOM residues into a consecutive sequence can corrupt BW offsets. Review engineered fusions, tags, unresolved helices, and long or ambiguous indels; leave uncertain regions unassigned or explicitly provisional.

Use evidence-backed `--anchors "3.50=131,6.50=288"`, `--tm5-end`, and `--tm6-start` only for local inference. These values are one-based positions in the full input sequence, not PDB residue IDs. Record the supporting evidence. Changing TM5/TM6 boundaries must not move 5.50 or 6.50.

## Deliver and check

The helper writes a complete TSV, readable HTML, and JSON containing the table and provenance. Show a concise preview with position, residue, `software_segment`, and `software_number`, including ECL2 rows around 45.50 and any requested loop/H8 regions, then link the full table. Include curated BW/GPCRdb columns when available and flag discrepancies. Do not replace a long table with an anchor-only result. Include one row per input sequence residue, even when BW is blank, and retain any extra unmapped coordinate residues. Keep unsupported positions blank rather than inventing labels.

Before delivery, verify row coverage, requested chain/record, original structure IDs, all software loop/H8 labels, ECL2 C45.50, insertion suffixes, and the seven helix .50 anchors where supported. Preserve the distinction between software inferences and curated labels. A curated table can include low-confidence software values marked `draft` in `software_number_status`; do not present them as validated. Low-confidence local or observed-only inputs do not export final tables by default. If an unresolved result is still useful, use `--draft`, keep `.draft` filenames and status, and explain what evidence would resolve it.

For changes to the helper or its engine, run `python scripts/test_annotation.py`. Do not alter the curated conventions to make a particular example pass.

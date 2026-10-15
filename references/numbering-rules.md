# Assignment rules and provenance

## Source implementation

Scientific source: **GPCRBWAnnotator v1.3**, `gpcr_bw_annotator.py`.

`scripts/bw_core.py` is the unchanged scientific prefix of `gpcr_bw_annotator.py`, ending before `annotation_tsv`. The original GUI and document writers are omitted. Original file SHA-256: `99889f10a0b1b3316fea6a41b7741674df6ae74e39d4b4cecbbb29dc910aced8`. Source version: 1.3.0; inspected 2026-10-15. The bundled engine is self-contained; the original installation and executable are not required. Do not run the desktop executable to annotate a sequence.

The source README says it aligns "locally"; the implementation actually performs **global Gotoh alignment**, BLOSUM62, gap opening 10 and extension 1. Follow the implementation.

## BW versus GPCRdb

For class-A TM helix H, BW is `H.(50 + p - a)`, where p is the target sequence position and a is that helix's .50 anchor position. Example: the residue immediately after 6.50 is 6.51. Counts are sequence-based, not author PDB residue arithmetic.

GPCRdb uses `x` for structure-based generic positions, correcting helix bulges and constrictions. A combined label such as `5.46x47` contains BW `5.46` and GPCRdb `5x47`. Never turn `5x47` into BW `5.47`; never derive a structural x-number from the local sequence-only engine. The curated API's `alternative_generic_numbers` entry with scheme `BW` is preferred; a verified class-A combined display label can supply its BW component.

BW is the native scheme for class A. Other classes use different native schemes; cross-class alternative BW entries do not make their native scheme BW.

## Bovine reference anchors and topology

Reference: canonical bovine rhodopsin, UniProt P02699, 348 aa. The source's topology is documented as curated from GPCRdb opsd_bovin on 2026-08-18.

| Label | Bovine residue | Segment |
|---|---|---|
| 1.50 | N55 | TM1 |
| 2.50 | D83 | TM2 |
| 3.50 | R135 | TM3 |
| 4.50 | W161 | TM4 |
| 45.50 | C187 | ECL2 |
| 5.50 | P215 | TM5 |
| 6.50 | P267 | TM6 |
| 7.50 | P303 | TM7 |

These are bovine template positions, never target positions copied without alignment. Local anchor candidates are scored within ±12 positions of the alignment projection, using identity, distance, and DRY-like, CWxP, and NPxxY context. ECL2 uses an independent cysteine selection based on relative position between 4.50 and 5.50, reducing cysteine-cluster errors. This is a heuristic, not proof of disulfide pairing.

Reference intervals: N-terminus 1–32; TM1 33–65; ICL1 66–69; TM2 70–101; ECL1 102–105; TM3 106–141; ICL2 142–148; TM4 149–174; ECL2 175–198; TM5 199–237; ICL3 238–239; TM6 240–278; ECL3 279–283; H8 310–322; TM7 284–309; C-terminus 323–348. Transfer intervals through alignment rather than applying reference residue numbers to the target.

## Full-loop and H8 convention

The complete v1.3 series is always returned in `software_number`, with `software_segment` identifying the region and `software_number_status` distinguishing an inference from a provisional value. Local inference also places loop/H8 labels in `loop_or_h8_number`, never `bw_number`:

| Loop | Segment code | Reference position receiving .50 |
|---|---|---|
| ICL1 | 12 | 68 |
| ECL1 | 23 | 103 |
| ICL2 | 34 | 142 |
| ECL2 | 45 | 187 |
| ICL3 | 56 | 239 |
| ECL3 | 67 | 281 |

Template-aligned positions retain their reference-relative labels. Target loop insertions receive alphabetic suffixes (`45.43a`, `45.43b`, …) so later homologous positions are not shifted. The complete ECL2 label set is adjusted around the independently detected C45.50. ICL3 residues moved by supported TM5/TM6 boundary overrides remain numbered relative to the projected 56.50 position. H8 starts at 8.47 and is alignment-relative in this engine; bovine F313 is 8.50.

GPCRdb has curated positions in some loops, but does not assign structural generic numbers to every flexible-loop residue. The v1.3 extension provides complete table labels, not universal structural equivalence. Curated mode still returns the complete v1.3 series in `software_number`; preserve the actual curated labels and blanks in their own columns. Never fill a curated structural x-number from a software loop label.

When a software anchor or boundary differs from the curated reference, retain both assignments and flag the affected positions in `summary.software_numbering.disagreements` and the row notes. This matters in cysteine-rich ECL2: the v1.3 cysteine-selection heuristic can select a different C45.50 from a receptor-specific curated assignment. Do not overwrite the curated anchor or silently merge two incompatible 45.xx series. A numbering discrepancy requires review before treating the sites as equivalent.

## Confidence and overrides

The source scores eight anchor identities (including C45.50), three motif contexts, anchor order, identity over the reference 7TM/H8 core, and TM hydrophobic fraction.

| Confidence | Minimum source criteria |
|---|---|
| High | 7/8 anchors; 2/3 motifs; ordered; core identity 28%; TM hydrophobicity 58% |
| Moderate | 6/8 anchors; 1/3 motifs; ordered; core identity 18%; TM hydrophobicity 52% |
| Low | Other results |

Source input range is 180–1,200 aa, one sequence per analysis. B, X, Z, J, U, and O are accepted; alignment gaps, position digits, whitespace, and terminal stop symbols are removed. Report ambiguous residues as supplied; do not replace them with a conserved anchor to increase confidence.

The helper adds draft-only treatment for observed-only structures and poor coordinate/reference mapping. A high local score cannot validate long fusions, divergent families, ambiguous loop alignments, partial receptors, or arbitrary non-GPCR proteins. Seek a curated alignment or structural evidence for these cases. Manual anchors are permitted by the original implementation and should carry independent evidence. Manual TM5 end must be at or after 5.50, TM6 start at or before 6.50, and TM5 end before TM6 start.

## Authoritative references

- [GPCRdb generic numbering](https://docs.gpcrdb.org/generic_numbering.html): BW/class-specific sequence schemes, x-numbering, bulges and constrictions.
- [GPCRdb web services](https://docs.gpcrdb.org/web_services.html): public protein and alignment API access.
- [Bovine rhodopsin](https://gpcrdb.org/protein/opsd_bovin/) and [UniProt P02699](https://www.uniprot.org/uniprotkb/P02699/entry).
- Ballesteros and Weinstein (1995), Methods in Neurosciences 25, 366–428; Isberg et al. (2015), Trends in Pharmacological Sciences 36, 22–31, as cited by the GPCRdb documentation.

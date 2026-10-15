# Input handling and output interpretation

## Structures

Use `--list` to inspect all chains, then select the receptor using PDB COMPND/entity descriptions and any verified UniProt/GPCRdb identity. Explicitly pass the selected author or label chain ID. The helper rejects ambiguous selection, including a chain identifier matching multiple records. It reads the first structural model and deduplicates alternate conformers at residue level. Author insertion codes, negative residue IDs, modified residues, and mmCIF label chain/sequence identifiers remain traceable.

For PDB, SEQRES defines the full construct sequence; align observed residues to it. For mmCIF, `_entity_poly_seq` defines the polymer sequence and `_atom_site.label_seq_id` directly maps observed residues without collapsing unresolved gaps. Residues outside the full sequence remain additional unassigned rows. Unknown protein residues are X; supported modified residues are mapped to the parent amino acid and noted in the table.

If only coordinates exist, fetch a matching full construct or canonical sequence and use:

```text
python scripts/annotate_gpcr.py --input receptor.pdb --chain A --reference-fasta full_sequence.fasta --out outputs/receptor_bw
```

Reference sequence positions and observed residue IDs are separate columns. Any residues added by the construct but absent from this full reference remain unmapped rows. Explain the chosen sequence basis. If a full sequence cannot be obtained, `--draft` permits a provisional observed-only table; compressed missing residues make its offsets unreliable.

Engineered fusion domains are not receptor loops. Prefer a matching curated receptor and leave fusion insertions unassigned. For a fusion with no curated match, isolate the verified GPCR sequence and retain an explicit map back to author residue IDs; otherwise deliver a provisional result. Do not present automatically filled local loop labels on a fusion protein as biological loop equivalence. Long ambiguous insertions, repeated sequences, and partial helices may require a curated alignment or structure-informed review even when overall identity is high.

## Remote and offline data

The helper accepts direct FASTA/PDB/mmCIF URLs, RCSB `/structure/ID`, UniProt `/uniprotkb/ACCESSION/entry`, and GPCRdb `/protein/entry_name/` pages. GPCRdb protein pages provide sequence; pass `--gpcrdb-entry` too to use curated residue numbers. Other HTML pages require explicit sequence extraction with browser/web tools or a linked data download. Authentication and transient site failures are not evidence about the protein.

Only public identifiers and URLs are fetched. Local inference does not upload sequences. If terminal network access is unavailable, use available browser/connector tools to retrieve the same source data and save it locally, or request the scoped network access supported by the environment. Never claim a fetch succeeded without checking the actual response.

For an offline curated run, place the unmodified JSON responses in a directory:

- `protein.json`: `https://gpcrdb.org/services/protein/ENTRY/`
- `residues.json`: `https://gpcrdb.org/services/residues/extended/ENTRY/`

Then add `--gpcrdb-entry ENTRY --gpcrdb-cache PATH`. Confirm entry name, species, and sequence. The helper validates endpoint agreement, native class-A scheme, and all residue positions. Keep retrieval dates and source URLs with the result if they are known.

Curated transfer uses a global alignment to map target sequence positions onto the selected reference. Identity below 90% or TM coverage below 80% requires an explicit draft. These are screening thresholds, not a validated homology detector. Verify the receptor identity separately. Substitutions retain the mapped reference label with a note. Target insertions have blank curated labels; the separate v1.3 `software_number` can supply an alignment-relative label, including alphabetic loop insertion suffixes. These are distinct claims. A closely matching engineered construct can retain curated receptor numbering across a deletion; its deleted positions are not residues in that input and do not receive synthetic rows.

Every run includes the full software numbering alongside the selected assignment method. With a curated reference, the software series is computed independently, so curated gaps cannot erase ECL2, other loop, or H8 labels. `summary.software_numbering` records its confidence, warnings, anchors, loop codes, and any disagreement with the curated numbering or topology. If software inference is low confidence, its values are marked `draft` even if the curated assignment is usable. If the software cannot analyze the input, its columns are blank and the summary records the reason; do not claim complete software numbering in that case.

## Output schema

Each output prefix yields `.tsv`, `.html`, and `.json`. Draft exports contain `.draft` before the extension. JSON contains `summary` and `residues`; HTML includes provenance and the full table. TSV is UTF-8 with BOM for desktop spreadsheet compatibility.

| Column | Meaning |
|---|---|
| sequence_position | One-based position in the declared or supplied sequence; blank for extra unmapped coordinate residues |
| residue | Amino acid in that sequence basis |
| chain / label_chain | Author chain ID and, for mmCIF, label chain ID |
| structure_residue | Original author residue number plus insertion code |
| observed_residue | Amino acid actually present in coordinates |
| observed | yes / no for structures; not applicable for sequence inputs |
| software_segment | Region assigned by the original v1.3 engine; can differ from curated topology |
| software_number | Complete original v1.3 BW / ECL1–3 / ICL1–3 / H8 label, including alphabetic loop insertion suffixes |
| software_number_status | inferred / draft / unassigned for the software series, independent of curated assignment status |
| segment | TM, loop, terminus, H8, or unassigned region |
| bw_number | Class-A sequence-based BW for TM1–TM7 only |
| gpcrdb_number | Curated structural x-number; blank for local inference |
| loop_or_h8_number | Curated non-TM display label or explicitly alignment-relative local label |
| assignment_source | Curated, transferred, v1.3 inference, or unassigned |
| status | assigned, draft, or unassigned |
| note | Source position, substitution, modification, or label interpretation |

Return all sequence residues, including termini and unresolved residues with blank BW. Include `software_number` in the table shown to the user; a BW-only or curated-only preview can hide the requested full-loop coverage. For structures, also retain extra observed residues that cannot be mapped. Do not use a PDB author's integer residue ID as an array index into a FASTA sequence. The summary records input/sequence hashes, numbering method, selected chain, confidence or mapping quality, and warnings. For local inference it also records all anchors and effective TM5/TM6 boundaries.

Exit code 0 means the selected operation completed; code 1 reports invalid input or a retrieval error; code 2 means uncertain inference was withheld from final export. For code 2, inspect the printed summary and resolve the cause or rerun with `--draft` when a provisional table serves the user's request.

# GPCR BW Numbering

A portable AI skill and standalone Python tool for residue-by-residue class-A GPCR numbering, based on **GPCRBWAnnotator v1.3**.

## What it provides

- Ballesteros–Weinstein (BW) numbering for TM1–TM7.
- Complete ECL2 **45.xx** numbering around the software's C45.50 anchor.
- ECL1/ECL3 labels **23.xx / 67.xx** and ICL1/ICL2/ICL3 labels **12.xx / 34.xx / 56.xx**.
- Helix-8 labels **8.xx** and loop insertion suffixes such as **45.43a**.
- Separate curated GPCRdb values, with disagreements flagged rather than silently merged.
- Full **TSV, HTML, and JSON** residue tables, preserving structure chain IDs, author residue numbers, insertion codes, and unresolved positions.

Accepts pasted sequences, FASTA, PDB/mmCIF (including gzip), PDB accessions, and protein links. Sequence and structure files are processed locally. Remote inputs and curated-reference mode retrieve public records; local sequences are not uploaded.

## Get the skill

Download the [ready-to-install skill ZIP](gpcr-bw-numbering.zip?raw=1) and extract its `gpcr-bw-numbering` folder.

Download this repository using **Code → Download ZIP**, extract it, and rename the extracted folder to `gpcr-bw-numbering`.

For Codex, place that folder in your personal skills directory, typically `~/.codex/skills/`. If `CODEX_HOME` is set, use its `skills/` directory. The installed folder must contain `SKILL.md`, `scripts/`, and `references/` together. Start a new chat and invoke:

```text
$gpcr-bw-numbering Number PDB 6OYA chain R and return a complete residue table, including ECL2, all other loops, and H8.
```

Other agents that support `SKILL.md` skills can use the same folder when their environment supports running Python. Follow the agent's own skill-import instructions; this repository does not require a Codex-specific API for computation.

## Standalone use

Requires **Python 3.10+**. No third-party Python packages are needed for annotation or the regression tests. Run commands from this repository's directory:

```bash
# Included public reference example
python scripts/annotate_gpcr.py --input examples/bovine-rhodopsin.fasta --out outputs/bovine

# Identify the receptor chain before selecting it
python scripts/annotate_gpcr.py --input receptor.pdb --list
python scripts/annotate_gpcr.py --input receptor.pdb --chain R --out outputs/receptor

# Fetch a public PDB entry (network access required)
python scripts/annotate_gpcr.py --pdb-id 6OYA --chain R --out outputs/6OYA_R

# Verified receptor-specific curated numbering plus the full software series
python scripts/annotate_gpcr.py --url https://www.uniprot.org/uniprotkb/P07550/entry --gpcrdb-entry adrb2_human --out outputs/adrb2
```

See `python scripts/annotate_gpcr.py --help` for pasted-sequence input, offline curated caches, full-sequence references for observed-only structures, manual anchors, and supported TM5/TM6 boundary overrides.

## Reading the table

`software_number` reproduces the complete v1.3 BW/loop/H8 series. `software_segment` identifies the region and `software_number_status` marks inferred, provisional, or unassigned values. Curated BW and GPCRdb structure-based `x` numbers remain in separate columns.

Full-loop and H8 labels from the software are **alignment-relative**. They do not establish structural equivalence in flexible loops. Curated GPCRdb and the software can disagree about an ECL2 anchor or segment boundary; those differences are retained in row notes and the JSON summary. Review them before comparing sites.

BW is the native sequence-based scheme for **class A**. Other GPCR classes require their appropriate schemes. A motif match or high heuristic confidence does not independently establish receptor identity. Engineered fusions, partial receptors, and ambiguous alignments need review. Observed-only structures require a matching full sequence to avoid compressing unresolved gaps into incorrect offsets. Low-confidence local results require `--draft` to export explicitly provisional files.

## Validation

```bash
python scripts/test_annotation.py
```

The suite covers known bovine rhodopsin and human β2-adrenoceptor anchors, complete loop/H8 output, insertion suffixes, supported ICL3 boundary changes, structure residue mapping, and draft handling. Development checks also verified all 348 original bovine sample labels and consistent mapping of all 326 observed receptor residues in the PDB/mmCIF versions of 6OYA.

## Scientific sources

- [GPCRdb generic numbering](https://docs.gpcrdb.org/generic_numbering.html)
- [GPCRdb web services](https://docs.gpcrdb.org/web_services.html)
- [Bovine rhodopsin, UniProt P02699](https://www.uniprot.org/uniprotkb/P02699/entry)
- Ballesteros and Weinstein (1995), Methods in Neurosciences 25, 366–428.
- Isberg et al. (2015), Trends in Pharmacological Sciences 36, 22–31.

Implementation provenance and detailed conventions are in [references/numbering-rules.md](references/numbering-rules.md). Input and output details are in [references/input-and-output.md](references/input-and-output.md).

## License

[MIT](LICENSE).

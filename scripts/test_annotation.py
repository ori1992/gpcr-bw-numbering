"""Regression tests for biological anchors, structural mapping, and draft handling."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import annotate_gpcr as app
import bw_core as core

ADRB2 = (
    'MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMW'
    'TFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYANETCCDFFTNQAY'
    'AIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVI'
    'QDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDN'
    'IDSQGRNCSTNDSLL')
THREE = {v: k for k, v in app.AA3.items() if k in app.AA3_STANDARD}


def args(**kwargs):
    values = dict(reference_fasta=None, gpcrdb_entry=None, gpcrdb_cache=None,
                  anchors='', tm5_end=None, tm6_start=None, draft=False)
    values.update(kwargs)
    return SimpleNamespace(**values)


def pdb_text(declared=True):
    lines = []
    sequence = core.BOVINE_SEQUENCE
    if declared:
        for start in range(0, len(sequence), 13):
            names = ' '.join(THREE[aa] for aa in sequence[start:start + 13])
            lines.append(f'SEQRES {start // 13 + 1:3d} R {len(sequence):4d}  {names}')
    for position, aa in enumerate(sequence, 1):
        if position in (80, 81):
            continue
        ins = 'A' if position == 101 else ' '
        lines.append(f'ATOM  {position:5d}  CA  {THREE[aa]:3s} R{position + 900:4d}{ins}   '
                     '   0.000   0.000   0.000  1.00 20.00           C')
    return '\n'.join(lines)


def cif_text():
    sequence = core.BOVINE_SEQUENCE
    lines = ['data_test', 'loop_', '_entity_poly_seq.entity_id', '_entity_poly_seq.num', '_entity_poly_seq.mon_id']
    lines.extend(f'1 {i} {THREE[aa]}' for i, aa in enumerate(sequence, 1))
    lines += ['#', 'loop_', '_struct_asym.id', '_struct_asym.entity_id', 'X 1', '#', 'loop_',
              '_atom_site.label_atom_id', '_atom_site.label_comp_id', '_atom_site.label_asym_id',
              '_atom_site.label_seq_id', '_atom_site.auth_asym_id', '_atom_site.auth_seq_id',
              '_atom_site.pdbx_PDB_ins_code', '_atom_site.pdbx_PDB_model_num']
    for i, aa in enumerate(sequence, 1):
        if i not in (80, 81):
            lines.append(f'CA {THREE[aa]} X {i} R {i + 900} ' + ('A' if i == 101 else '?') + ' 1')
    lines.append('#')
    return '\n'.join(lines)


class AnnotationTests(unittest.TestCase):
    def test_bovine_anchors_and_complete_loop_labels(self):
        protein = app.parse_fasta('>bovine\n' + core.BOVINE_SEQUENCE)[0]
        rows, summary, doubtful = app.build_table(protein, args())
        self.assertFalse(doubtful)
        self.assertEqual(summary['confidence'], 'High')
        self.assertEqual(len(rows), 348)
        for label, position, aa, segment in core.ANCHORS:
            row = rows[position - 1]
            self.assertEqual(row['residue'], aa)
            self.assertEqual(row['segment'], segment)
            self.assertEqual(row['bw_number'] if segment.startswith('TM') else row['loop_or_h8_number'], label)
        for segment in core.LOOP_CODES:
            self.assertTrue(all(r['loop_or_h8_number'] and not r['bw_number'] for r in rows if r['segment'] == segment))
        original = core.analyze_sequence(core.SequenceInput('bovine', core.BOVINE_SEQUENCE))
        self.assertEqual([(r['software_segment'], r['software_number']) for r in rows],
                         [(r.segment, r.generic_number) for r in original.annotations])
        self.assertEqual(rows[312]['loop_or_h8_number'], '8.50')
        self.assertFalse(any(r['gpcrdb_number'] for r in rows))

    def test_non_opsin_reference_anchors(self):
        result = core.analyze_sequence(core.SequenceInput('ADRB2_HUMAN', ADRB2))
        self.assertEqual({a.label: a.position for a in result.anchors},
                         {'1.50': 51, '2.50': 79, '3.50': 131, '4.50': 158,
                          '45.50': 184, '5.50': 211, '6.50': 288, '7.50': 323})
        self.assertIn(result.confidence, ('High', 'Moderate'))

    def test_loop_insertions_keep_homologous_anchor(self):
        sequence = core.BOVINE_SEQUENCE[:180] + 'AA' + core.BOVINE_SEQUENCE[180:238] + 'GGG' + core.BOVINE_SEQUENCE[238:]
        result = core.analyze_sequence(core.SequenceInput('insertions', sequence))
        self.assertEqual(next(r.position for r in result.annotations if r.generic_number == '45.50'), 189)
        self.assertGreaterEqual(sum(r.generic_number.startswith('45.') and r.generic_number[-1:].isalpha()
                                    for r in result.annotations), 2)

    def test_manual_boundaries_preserve_anchors(self):
        result = core.analyze_sequence(core.SequenceInput('bovine', core.BOVINE_SEQUENCE),
                                       manual_tm5_end=230, manual_tm6_start=241)
        self.assertEqual({a.label: a.position for a in result.anchors}, {l: p for l, p, *_ in core.ANCHORS})
        self.assertEqual(result.annotations[230].segment, 'ICL3')
        self.assertEqual(result.annotations[240].generic_number, '6.24')
        with self.assertRaises(ValueError):
            core.analyze_sequence(core.SequenceInput('bovine', core.BOVINE_SEQUENCE),
                                  manual_tm5_end=250, manual_tm6_start=240)

    def check_structure(self, protein):
        rows, summary, doubtful = app.build_table(protein, args())
        self.assertFalse(doubtful)
        self.assertEqual(len(rows), 348)
        self.assertEqual(rows[79]['observed'], 'no')
        self.assertEqual(rows[80]['observed'], 'no')
        self.assertEqual(rows[82]['structure_residue'], '983')
        self.assertEqual(rows[82]['bw_number'], '2.50')
        self.assertEqual(rows[100]['structure_residue'], '1001A')
        self.assertEqual(summary['observed_residues'], 346)

    def test_pdb_missing_positions_and_insertion_codes(self):
        self.check_structure(app.parse_pdb(pdb_text())[0])

    def test_mmcif_direct_polymer_mapping(self):
        protein = app.parse_mmcif(cif_text())[0]
        self.assertEqual((protein.chain, protein.label_chain), ('R', 'X'))
        self.check_structure(protein)

    def test_multichain_selection_is_explicit(self):
        records = [app.Protein('first', core.BOVINE_SEQUENCE, 'A'), app.Protein('second', core.BOVINE_SEQUENCE, 'B')]
        with self.assertRaisesRegex(ValueError, 'Select one record'):
            app.choose(records)
        self.assertEqual(app.choose(records, 'B').chain, 'B')

    def test_fasta_records_and_invalid_html(self):
        records = app.parse_fasta('>one\n' + ADRB2 + '\n>two\n' + core.BOVINE_SEQUENCE)
        self.assertEqual(app.choose(records, record='two').sequence, core.BOVINE_SEQUENCE)
        with self.assertRaises(ValueError):
            app.parse_fasta('<html>protein sequence</html>')

    def test_observed_only_is_draft_and_reference_repairs_missing_positions(self):
        protein = app.parse_pdb(pdb_text(False))[0]
        _, _, doubtful = app.build_table(protein, args())
        self.assertTrue(doubtful)
        with tempfile.TemporaryDirectory() as tmp:
            full = Path(tmp) / 'full.fasta'
            full.write_text(core.BOVINE_SEQUENCE)
            rows, _, doubtful = app.build_table(protein, args(reference_fasta=full))
            self.assertFalse(doubtful)
            self.assertEqual(rows[82]['structure_residue'], '983')
            self.assertEqual(rows[82]['bw_number'], '2.50')

    def test_curated_bw_and_structural_number_are_distinct(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            info = dict(entry_name='opsd_bovin', sequence=core.BOVINE_SEQUENCE,
                        receptor_class='Class A (Rhodopsin)', residue_numbering_scheme='GPCRdb(A)')
            result = core.analyze_sequence(core.SequenceInput('bovine', core.BOVINE_SEQUENCE))
            residues = [dict(sequence_number=r.position, amino_acid=r.residue, protein_segment=r.segment,
                             display_generic_number=None, alternative_generic_numbers=[]) for r in result.annotations]
            residues[210].update(display_generic_number='5.46x47',
                                 alternative_generic_numbers=[dict(scheme='BW', label='5.46')])
            (folder / 'protein.json').write_text(json.dumps(info))
            (folder / 'residues.json').write_text(json.dumps(residues))
            rows, _, doubtful = app.build_table(app.Protein('bovine', core.BOVINE_SEQUENCE),
                                                args(gpcrdb_entry='opsd_bovin', gpcrdb_cache=folder))
            self.assertFalse(doubtful)
            self.assertEqual((rows[210]['bw_number'], rows[210]['gpcrdb_number']), ('5.46', '5x47'))
            for segment in (*core.LOOP_CODES, 'Helix 8'):
                selected = [r for r in rows if r['software_segment'] == segment]
                self.assertTrue(selected, segment)
                self.assertTrue(all(r['software_number'] and r['software_number_status'] == 'inferred' for r in selected), segment)
                self.assertTrue(all(not r['gpcrdb_number'] and not r['loop_or_h8_number'] for r in selected), segment)
            self.assertEqual(rows[186]['software_number'], '45.50')
            self.assertEqual(rows[312]['software_number'], '8.50')
            self.assertEqual(rows[210]['software_number'], '5.46')
            # Deliberately conflicting curated fixture: neither numbering may overwrite the other.
            residues[186]['display_generic_number'] = '45.49x49'
            (folder / 'residues.json').write_text(json.dumps(residues))
            rows, summary, _ = app.build_table(app.Protein('bovine', core.BOVINE_SEQUENCE),
                                               args(gpcrdb_entry='opsd_bovin', gpcrdb_cache=folder))
            self.assertEqual((rows[186]['software_number'], rows[186]['loop_or_h8_number']), ('45.50', '45.49x49'))
            self.assertTrue(any(r['sequence_position'] == 187 for r in summary['software_numbering']['disagreements']))
            info['receptor_class'] = 'Class B1 (Secretin)'
            (folder / 'protein.json').write_text(json.dumps(info))
            with self.assertRaisesRegex(ValueError, 'not class A'):
                app.curated_annotations(app.Protein('bovine', core.BOVINE_SEQUENCE), 'opsd_bovin', folder)

    def test_software_insertions_and_manual_icl3_are_exported(self):
        sequence = core.BOVINE_SEQUENCE[:180] + 'AA' + core.BOVINE_SEQUENCE[180:238] + 'GGG' + core.BOVINE_SEQUENCE[238:]
        result = core.analyze_sequence(core.SequenceInput('insertions', sequence))
        rows, _, _ = app.build_table(app.Protein('insertions', sequence), args())
        self.assertEqual([r['software_number'] for r in rows], [r.generic_number for r in result.annotations])
        self.assertTrue(any(r['software_number'].startswith('45.') and r['software_number'][-1:].isalpha() for r in rows))
        rows, _, _ = app.build_table(app.Protein('bovine', core.BOVINE_SEQUENCE), args(tm5_end=230, tm6_start=241))
        self.assertTrue(all(r['software_segment'] == 'ICL3' and r['software_number'].startswith('56.') for r in rows[230:240]))
        self.assertEqual(rows[240]['software_number'], '6.24')

    def test_cli_exports_complete_table_and_blocks_final_low_confidence(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            prefix = str(Path(tmp) / 'result')
            self.assertEqual(app.main(['--sequence', core.BOVINE_SEQUENCE, '--out', prefix]), 0)
            document = json.loads(Path(prefix + '.json').read_text(encoding='utf-8'))
            self.assertEqual(len(document['residues']), 348)
            self.assertEqual(document['residues'][186]['software_number'], '45.50')
            self.assertEqual(len(Path(prefix + '.tsv').read_text(encoding='utf-8-sig').splitlines()), 349)
            low = str(Path(tmp) / 'low')
            self.assertEqual(app.main(['--sequence', 'A' * 348, '--out', low]), 2)
            self.assertFalse(Path(low + '.tsv').exists())
            self.assertEqual(app.main(['--sequence', 'A' * 348, '--out', low, '--draft']), 0)
            self.assertTrue(Path(low + '.draft.tsv').exists())
            provisional = json.loads(Path(low + '.draft.json').read_text(encoding='utf-8'))
            self.assertTrue(all(r['software_number_status'] == 'draft' for r in provisional['residues'] if r['software_number']))

    def test_known_website_links(self):
        self.assertEqual(app.resolve_url('https://www.rcsb.org/structure/6OYA'), 'https://files.rcsb.org/download/6OYA.cif')
        self.assertEqual(app.resolve_url('https://www.uniprot.org/uniprotkb/P07550/entry'), 'https://rest.uniprot.org/uniprotkb/P07550.fasta')
        self.assertEqual(app.resolve_url('https://gpcrdb.org/protein/adrb2_human/'), 'https://gpcrdb.org/services/protein/adrb2_human/')


if __name__ == '__main__':
    unittest.main(verbosity=2)

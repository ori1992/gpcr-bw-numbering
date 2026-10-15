#!/usr/bin/env python3
"""Create traceable class-A GPCR residue tables. Python 3.10+, standard library only."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import io
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import bw_core as core

AA3 = dict(zip(
    'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASX GLX UNK SEC PYL'.split(),
    'ARNDCQEGHILKMFPSTWYVBZXUO'))
AA3.update(MSE='M', SEP='S', TPO='T', PTR='Y', HYP='P', CSO='C')
FIELDS = ['sequence_position', 'residue', 'chain', 'label_chain', 'structure_residue',
          'observed_residue', 'observed', 'software_segment', 'software_number',
          'software_number_status', 'segment', 'bw_number', 'gpcrdb_number',
          'loop_or_h8_number', 'assignment_source', 'status', 'note']


@dataclass
class Observed:
    residue: str
    auth_id: str
    chain: str
    label_chain: str = ''
    label_position: int | None = None
    original_name: str = ''


@dataclass
class Protein:
    name: str
    sequence: str
    chain: str = ''
    label_chain: str = ''
    sequence_source: str = 'sequence'
    observed: list[Observed] = field(default_factory=list)


def alignment_map(reference, target):
    """Map target positions onto reference positions; never close a reference gap."""
    left, right, _ = core.pairwise_align(reference, target)
    mapping, p, q, paired, matches = {}, 0, 0, 0, 0
    for a, b in zip(left, right):
        p += a != '-'
        q += b != '-'
        if b != '-':
            mapping[q] = p if a != '-' else None
        if a != '-' and b != '-':
            paired += 1
            matches += a == b
    return mapping, matches / paired if paired else 0.0


def parse_fasta(text):
    if text.lstrip().startswith('<'):
        raise ValueError('HTML is not a sequence. Extract a FASTA or sequence field from the linked page.')
    records, header, lines = [], None, []
    for line in text.strip().splitlines():
        if line.lstrip().startswith('>'):
            if lines:
                parsed = core.parse_sequence_text('\n'.join(lines), header or 'GPCR')
                records.append(Protein(parsed.name, parsed.sequence))
            header, lines = line.strip()[1:].strip(), []
        elif line.strip():
            lines.append(line)
    if lines:
        parsed = core.parse_sequence_text('\n'.join(lines), header or 'GPCR')
        records.append(Protein(parsed.name, parsed.sequence))
    if not records:
        raise ValueError('No sequence records found.')
    return records


def parse_pdb(text):
    declared, atoms, seen = {}, {}, set()
    model_seen = False
    for line in text.splitlines():
        kind = line[:6].strip()
        if kind == 'SEQRES':
            declared.setdefault(line[11:12].strip(), []).extend(line[19:70].split())
        if kind == 'MODEL':
            if model_seen:
                break
            model_seen = True
        if kind == 'ENDMDL':
            break
        if kind not in ('ATOM', 'HETATM') or line[12:16].strip() != 'CA':
            continue
        name = line[17:20].strip()
        if kind == 'HETATM' and name not in AA3:
            continue
        chain, residue_id = line[21:22].strip(), line[22:26].strip() + line[26:27].strip()
        key = (chain, residue_id)
        if key in seen:
            continue
        seen.add(key)
        atoms.setdefault(chain, []).append(Observed(AA3.get(name, 'X'), residue_id, chain, original_name=name))
    records = []
    for chain in dict.fromkeys([*declared, *atoms]):
        observed = atoms.get(chain, [])
        sequence = ''.join(AA3.get(aa, 'X') for aa in declared.get(chain, []))
        source = 'PDB SEQRES' if sequence else 'observed coordinates only'
        sequence = sequence or ''.join(item.residue for item in observed)
        if sequence:
            records.append(Protein('chain ' + (chain or '<blank>'), sequence, chain,
                                   sequence_source=source, observed=observed))
    if not records:
        raise ValueError('No protein chains found in PDB.')
    return records


def cif_tokens(text):
    """CIF 1.1 tokenization: quotes, comments, and column-one semicolon fields."""
    lines = iter(text.splitlines())
    for line in lines:
        if line.startswith(';'):
            block = [line[1:]]
            for following in lines:
                if following.startswith(';'):
                    break
                block.append(following)
            else:
                raise ValueError('Unterminated mmCIF multiline value.')
            yield '\n'.join(block), True
            continue
        i = 0
        while i < len(line):
            while i < len(line) and line[i].isspace():
                i += 1
            if i == len(line) or line[i] == '#':
                break
            start, quote = i, line[i]
            if quote in ('"', "'"):
                i += 1
                start = i
                while i < len(line):
                    if line[i] == quote and (i + 1 == len(line) or line[i + 1].isspace()):
                        break
                    i += 1
                if i == len(line):
                    raise ValueError('Unterminated mmCIF quoted value.')
                yield line[start:i], True
                i += 1
            else:
                while i < len(line) and not line[i].isspace():
                    i += 1
                yield line[start:i], False


def cif_categories(text):
    tokens, categories, i = list(cif_tokens(text)), {}, 0
    def control(token):
        value, quoted = token
        return not quoted and (value.startswith(('_', 'data_', 'save_')) or value in ('loop_', 'stop_', 'global_'))
    while i < len(tokens):
        value, quoted = tokens[i]
        if value == 'loop_' and not quoted:
            i += 1
            keys = []
            while i < len(tokens) and not tokens[i][1] and tokens[i][0].startswith('_'):
                keys.append(tokens[i][0])
                i += 1
            if not keys:
                raise ValueError('mmCIF loop without column names.')
            values = []
            while i < len(tokens) and not (len(values) % len(keys) == 0 and control(tokens[i])):
                values.append(tokens[i][0])
                i += 1
            if len(values) % len(keys):
                raise ValueError('Incomplete mmCIF loop row.')
            category = keys[0].split('.')[0]
            if category in ('_entity_poly_seq', '_struct_asym', '_atom_site'):
                categories.setdefault(category, []).extend(
                    dict(zip(keys, values[j:j + len(keys)])) for j in range(0, len(values), len(keys)))
        elif not quoted and value.startswith('_'):
            if i + 1 >= len(tokens):
                raise ValueError('Missing mmCIF item value.')
            category = value.split('.')[0]
            rows = categories.setdefault(category, [{}])
            rows[0][value] = tokens[i + 1][0]
            i += 2
        else:
            i += 1
    return categories


def parse_mmcif(text):
    cats = cif_categories(text)
    entities = {}
    for row in cats.get('_entity_poly_seq', []):
        entity = row['_entity_poly_seq.entity_id']
        position = int(row['_entity_poly_seq.num'])
        entities.setdefault(entity, {}).setdefault(position, AA3.get(row['_entity_poly_seq.mon_id'], 'X'))
    chains = {row['_struct_asym.id']: row['_struct_asym.entity_id'] for row in cats.get('_struct_asym', [])}
    atoms, auth_chains, seen, first_model = {}, {}, set(), None
    for row in cats.get('_atom_site', []):
        atom = row.get('_atom_site.label_atom_id', row.get('_atom_site.auth_atom_id', ''))
        name = row.get('_atom_site.label_comp_id', '')
        label_position = row.get('_atom_site.label_seq_id', '.')
        if atom != 'CA' or label_position in ('.', '?'):
            continue
        model = row.get('_atom_site.pdbx_PDB_model_num', '1')
        first_model = first_model or model
        if model != first_model:
            continue
        label_chain = row['_atom_site.label_asym_id']
        chain = row.get('_atom_site.auth_asym_id', label_chain)
        auth_chains[label_chain] = chain
        ins = row.get('_atom_site.pdbx_PDB_ins_code', '?')
        residue_id = row.get('_atom_site.auth_seq_id', label_position) + ('' if ins in ('.', '?') else ins)
        key = (label_chain, int(label_position))
        if key in seen:
            continue
        seen.add(key)
        atoms.setdefault(label_chain, []).append(Observed(
            AA3.get(name, 'X'), residue_id, chain, label_chain, int(label_position), name))
    records = []
    for label_chain in dict.fromkeys([*chains, *atoms]):
        positions = entities.get(chains.get(label_chain), {})
        observed = atoms.get(label_chain, [])
        # Keep protein polymers only (nucleic acids contain no recognized AA3 monomers).
        if not observed and (not positions or all(aa == 'X' for aa in positions.values())):
            continue
        sequence = ''.join(positions.get(i, 'X') for i in range(1, max(positions, default=0) + 1))
        source = 'mmCIF entity_poly_seq' if sequence else 'observed coordinates only'
        sequence = sequence or ''.join(item.residue for item in observed)
        chain = auth_chains.get(label_chain, label_chain)
        records.append(Protein('chain ' + chain, sequence, chain, label_chain, source, observed))
    if not records:
        raise ValueError('No protein polymer chains found in mmCIF.')
    return records


def fetch(url):
    if urlparse(url).scheme not in ('https', 'http'):
        raise ValueError('Use an HTTP(S) data URL or a local file.')
    request = Request(url, headers={'User-Agent': 'GPCR-BW-Skill/1.0'})
    with urlopen(request, timeout=30) as response:
        data = response.read(30_000_001)
    if len(data) > 30_000_000:
        raise ValueError('Remote input exceeds 30 MB; download a selected protein chain locally.')
    return data


def resolve_url(url):
    parsed = urlparse(url)
    if parsed.hostname in ('rcsb.org', 'www.rcsb.org'):
        match = re.search(r'/(?:structure|explore/explore.do)/([A-Za-z0-9]+)', parsed.path)
        if match:
            return 'https://files.rcsb.org/download/' + match.group(1).upper() + '.cif'
    if parsed.hostname in ('uniprot.org', 'www.uniprot.org'):
        match = re.search(r'/uniprot(?:kb)?/([A-Za-z0-9-]+)', parsed.path)
        if match:
            return 'https://rest.uniprot.org/uniprotkb/' + match.group(1) + '.fasta'
    if parsed.hostname in ('gpcrdb.org', 'www.gpcrdb.org'):
        match = re.fullmatch(r'/protein/([a-z0-9_]+)/?', parsed.path)
        if match:
            return 'https://gpcrdb.org/services/protein/' + match.group(1) + '/'
    return url


def decode_records(data):
    if data[:2] == b'\x1f\x8b':
        data = gzip.decompress(data)
    text = data.decode('utf-8-sig')
    if text.lstrip().startswith('data_'):
        return parse_mmcif(text)
    if any(line.startswith(('ATOM  ', 'SEQRES', 'HEADER', 'HETATM')) for line in text.splitlines()):
        return parse_pdb(text)
    if text.lstrip().startswith('{'):
        info = json.loads(text)
        if info.get('receptor_class') and not info['receptor_class'].startswith('Class A'):
            raise ValueError('The linked receptor is not class A; BW is not its native numbering scheme.')
        if not isinstance(info.get('sequence'), str):
            raise ValueError('JSON must contain an explicit protein sequence field.')
        return parse_fasta('>' + info.get('entry_name', 'linked_protein') + '\n' + info['sequence'])
    return parse_fasta(text)


def choose(records, chain=None, record=None):
    if chain is not None:
        selected = [r for r in records if r.chain == chain or r.label_chain == chain]
    elif record:
        selected = [r for r in records if r.name == record or r.name.split()[0] == record]
    else:
        selected = records
    if len(selected) != 1:
        options = [{'record': r.name, 'chain': r.chain, 'label_chain': r.label_chain,
                    'length': len(r.sequence)} for r in records]
        raise ValueError('Select one record using --chain or --record. Available: ' + json.dumps(options))
    return selected[0]


def observation_map(protein):
    if not protein.observed:
        return {}, [], []
    warnings = []
    if protein.sequence_source == 'mmCIF entity_poly_seq':
        mapping = {i: item.label_position for i, item in enumerate(protein.observed, 1)}
    elif protein.sequence_source == 'observed coordinates only':
        mapping = {i: i for i in range(1, len(protein.observed) + 1)}
    else:
        mapping, identity = alignment_map(protein.sequence, ''.join(item.residue for item in protein.observed))
        if identity < .90:
            warnings.append('Coordinate-to-sequence alignment identity is below 90%; inspect mapping.')
    assigned, extras = {}, []
    for i, observed in enumerate(protein.observed, 1):
        position = mapping.get(i)
        if position is None or not 1 <= position <= len(protein.sequence):
            extras.append(observed)
        elif position in assigned:
            raise ValueError('Multiple observed residues map to the same sequence position.')
        else:
            assigned[position] = observed
    return assigned, extras, warnings


def curated_annotations(protein, entry, cached=None):
    if not re.fullmatch(r'[a-z0-9]+_[a-z0-9]+', entry):
        raise ValueError('Use a GPCRdb entry name such as adrb2_human.')
    base = 'https://gpcrdb.org/services/'
    if cached:
        info = json.loads((cached / 'protein.json').read_text(encoding='utf-8'))
        residues = json.loads((cached / 'residues.json').read_text(encoding='utf-8'))
    else:
        info = json.loads(fetch(base + 'protein/' + entry + '/'))
        residues = json.loads(fetch(base + 'residues/extended/' + entry + '/'))
    if info.get('entry_name') != entry:
        raise ValueError('Curated protein cache does not match the requested entry.')
    if info.get('residue_numbering_scheme') != 'GPCRdb(A)' or not info.get('receptor_class', '').startswith('Class A'):
        raise ValueError('This entry is not class A. BW is a class-A scheme; report the appropriate class-specific scheme.')
    reference = info['sequence']
    lookup = {int(r['sequence_number']): r for r in residues}
    if len(lookup) != len(residues) or sorted(lookup) != list(range(1, len(reference) + 1)):
        raise ValueError('Curated residue records are incomplete or contain duplicate sequence positions.')
    if ''.join(lookup[i]['amino_acid'] for i in range(1, len(reference) + 1)) != reference:
        raise ValueError('Curated sequence and residue endpoint disagree.')
    mapping, identity = alignment_map(reference, protein.sequence)
    numbered = {p for p, r in lookup.items() if re.fullmatch(r'TM[1-7]', r['protein_segment'])}
    coverage = len(numbered.intersection(p for p in mapping.values() if p)) / max(1, len(numbered))
    doubtful = identity < .90 or coverage < .80
    annotations, warnings = [], []
    if doubtful:
        warnings.append('Curated reference mapping is below 90% identity or 80% TM coverage; draft only.')
    for position in range(1, len(protein.sequence) + 1):
        item = lookup.get(mapping.get(position))
        if item is None:
            annotations.append(dict(segment='Unassigned insertion', bw_number='', gpcrdb_number='',
                                    loop_or_h8_number='', assignment_source='unassigned',
                                    note='Target insertion has no curated reference equivalent.'))
            continue
        segment = item['protein_segment']
        display = item.get('display_generic_number') or ''
        alternatives = {r['scheme']: r['label'] for r in item.get('alternative_generic_numbers', [])}
        bw = alternatives.get('BW', '')
        if not bw and re.fullmatch(r'[1-7]\.\d+x\d+', display):
            bw = display.split('x')[0]
        gpcrdb = (display.split('.')[0] + 'x' + display.split('x')[1]) if 'x' in display else ''
        note = 'GPCRdb ' + entry + ', reference position ' + str(mapping[position])
        if protein.sequence[position - 1] != item['amino_acid']:
            note += '; aligned substitution ' + item['amino_acid'] + '→' + protein.sequence[position - 1]
        annotations.append(dict(segment=segment, bw_number=bw if re.fullmatch(r'TM[1-7]', segment) else '',
                                gpcrdb_number=gpcrdb, loop_or_h8_number=display if not re.fullmatch(r'TM[1-7]', segment) else '',
                                assignment_source='curated GPCRdb' if reference == protein.sequence else 'curated-reference transfer',
                                note=note))
    return annotations, dict(method='GPCRdb curated reference', reference=entry,
                             reference_identity=identity, reference_tm_coverage=coverage,
                             confidence='Low' if doubtful else 'Curated reference mapping',
                             warnings=warnings, endpoints=[base + 'protein/' + entry + '/',
                                                          base + 'residues/extended/' + entry + '/']), doubtful


def local_annotations(protein, args):
    parsed = core.parse_sequence_text(protein.sequence, protein.name)
    result = core.analyze_sequence(parsed, args.anchors, manual_tm5_end=args.tm5_end,
                                   manual_tm6_start=args.tm6_start)
    annotations = []
    for item in result.annotations:
        tm = bool(re.fullmatch(r'TM[1-7]', item.segment))
        annotations.append(dict(segment=item.segment, bw_number=item.generic_number if tm else '',
                                gpcrdb_number='', loop_or_h8_number=item.generic_number if not tm else '',
                                assignment_source='GPCRBWAnnotator 1.3 inference', note=item.annotation))
    return annotations, dict(method='GPCRBWAnnotator 1.3 rhodopsin alignment', reference='P02699',
                             confidence=result.confidence, confidence_detail=result.confidence_detail,
                             anchors=[asdict(a) for a in result.anchors],
                             topology_boundaries=asdict(result.topology_boundaries),
                             warnings=list(result.warnings)), result.confidence == 'Low'


def include_software_numbering(protein, annotations, summary, curated=False):
    """Always expose the complete v1.3 table without overwriting curated labels."""
    if curated:
        try:
            parsed = core.parse_sequence_text(protein.sequence, protein.name)
            result = core.analyze_sequence(parsed)
        except ValueError as exc:
            for annotation in annotations:
                annotation.update(software_segment='', software_number='', software_number_status='unassigned')
            summary['software_numbering'] = dict(engine='GPCRBWAnnotator 1.3', confidence='Unavailable',
                                                  warnings=[str(exc)], disagreements=[])
            summary['warnings'].append('Full software numbering unavailable: ' + str(exc))
            return
        software = [(item.segment, item.generic_number) for item in result.annotations]
        confidence = result.confidence
        detail = result.confidence_detail
        warnings = list(result.warnings)
        anchors = [asdict(item) for item in result.anchors]
    else:
        software = [(item['segment'], item['bw_number'] or item['loop_or_h8_number']) for item in annotations]
        confidence = summary['confidence']
        detail = summary['confidence_detail']
        warnings = list(summary['warnings'])
        anchors = summary['anchors']
    disagreements = []
    for position, (annotation, (segment, number)) in enumerate(zip(annotations, software), 1):
        annotation.update(software_segment=segment, software_number=number,
                          software_number_status=('draft' if confidence == 'Low' else 'inferred') if number else 'unassigned')
        if not number:
            continue
        if segment in core.LOOP_CODES or segment == 'Helix 8':
            annotation['note'] += '; full v1.3 ' + segment + ' label ' + number + ' is alignment-relative'
        if not curated:
            continue
        curated_segment = {'H8': 'Helix 8', 'N-term': 'N-terminus', 'C-term': 'C-terminus'}.get(annotation['segment'], annotation['segment'])
        curated_number = (annotation['bw_number'] or annotation['loop_or_h8_number']).split('x')[0]
        if curated_segment != segment or (curated_number and curated_number != number):
            disagreements.append(dict(sequence_position=position, residue=protein.sequence[position - 1],
                                      software_segment=segment, software_number=number,
                                      curated_segment=curated_segment, curated_number=curated_number))
            annotation['note'] += '; v1.3 numbering/topology differs from curated reference; review before comparing sites'
    summary['software_numbering'] = dict(engine='GPCRBWAnnotator 1.3', confidence=confidence,
                                          confidence_detail=detail, warnings=warnings, anchors=anchors,
                                          loop_codes=core.LOOP_CODES, h8_start='8.47', disagreements=disagreements)
    if disagreements:
        summary['warnings'].append(str(len(disagreements)) + ' software numbering/topology differences from curated reference; see software_numbering.disagreements.')


def build_table(protein, args):
    if args.reference_fasta:
        if protein.sequence_source != 'observed coordinates only':
            raise ValueError('--reference-fasta is for structures without a declared polymer sequence.')
        reference = choose(parse_fasta(args.reference_fasta.read_text(encoding='utf-8-sig')))
        protein.sequence, protein.sequence_source = reference.sequence, 'supplied full reference sequence'
    if args.gpcrdb_entry:
        if args.anchors or args.tm5_end or args.tm6_start:
            raise ValueError('Manual anchors/boundaries apply to local inference, not curated reference transfer.')
        annotations, summary, doubtful = curated_annotations(protein, args.gpcrdb_entry, args.gpcrdb_cache)
    else:
        annotations, summary, doubtful = local_annotations(protein, args)
    include_software_numbering(protein, annotations, summary, curated=bool(args.gpcrdb_entry))
    observed, extras, warnings = observation_map(protein)
    summary['warnings'].extend(warnings)
    doubtful |= bool(warnings)
    if protein.sequence_source == 'observed coordinates only':
        summary['warnings'].append('No full polymer sequence: missing residues may compress BW offsets. Supply --reference-fasta. Draft only.')
        doubtful = True
    status = 'draft' if doubtful or args.draft else 'assigned'
    rows = []
    for position, (residue, annotation) in enumerate(zip(protein.sequence, annotations), 1):
        atom = observed.get(position)
        row = dict(sequence_position=position, residue=residue, chain=protein.chain, label_chain=protein.label_chain,
                   structure_residue=atom.auth_id if atom else '', observed_residue=atom.residue if atom else '',
                   observed=('yes' if atom else 'no') if protein.sequence_source != 'sequence' else 'not applicable',
                   status='unassigned' if annotation['assignment_source'] == 'unassigned' else status, **annotation)
        if status == 'draft' and row['software_number']:
            row['software_number_status'] = 'draft'
        if atom and atom.residue != residue:
            row['note'] += '; coordinate residue differs from sequence'
        if atom and atom.original_name not in AA3_STANDARD:
            row['note'] += '; modified residue ' + atom.original_name
        rows.append(row)
    for atom in extras:
        rows.append(dict(sequence_position='', residue=atom.residue, chain=atom.chain,
                         label_chain=atom.label_chain, structure_residue=atom.auth_id, observed_residue=atom.residue,
                         observed='yes', segment='Unmapped coordinate residue', bw_number='', gpcrdb_number='',
                         software_segment='', software_number='', software_number_status='unassigned',
                         loop_or_h8_number='', assignment_source='unassigned', status='unassigned',
                         note='No full-sequence equivalent; residue retained without invented numbering.'))
    summary.update(sequence_name=protein.name, chain=protein.chain, label_chain=protein.label_chain,
                   sequence_source=protein.sequence_source, sequence_length=len(protein.sequence),
                   sequence_sha256=hashlib.sha256(protein.sequence.encode()).hexdigest(),
                   observed_residues=len(protein.observed), unmapped_observed=len(extras), status=status,
                   numbering_note='software_number reproduces the complete v1.3 BW/loop/H8 labels and insertion suffixes; loop/H8 labels are alignment-relative. Curated BW and GPCRdb x-numbers remain separate.')
    return rows, summary, doubtful


AA3_STANDARD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())


def export(rows, summary, prefix, draft):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    stem = str(prefix) + ('.draft' if draft else '')
    paths = [Path(stem + ext) for ext in ('.tsv', '.html', '.json')]
    with paths[0].open('w', encoding='utf-8-sig', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    esc = lambda value: html.escape(str(value))
    head = ''.join('<th>' + esc(key) + '</th>' for key in FIELDS)
    body = ''.join('<tr>' + ''.join('<td>' + esc(row.get(key, '')) + '</td>' for key in FIELDS) + '</tr>' for row in rows)
    content = ('<!doctype html><html lang="en"><meta charset="utf-8"><title>GPCR BW residue table</title>'
               '<style>body{font:14px system-ui;margin:24px}table{border-collapse:collapse}th,td{border:1px solid #ccc;'
               'padding:5px;text-align:left}th{position:sticky;top:0;background:#e7eef5}tr:nth-child(even){background:#f4f6f8}'
               'pre{white-space:pre-wrap}</style><h1>GPCR BW residue table' + (' — DRAFT' if draft else '') + '</h1>'
               '<pre>' + esc(json.dumps(summary, indent=2, ensure_ascii=False)) + '</pre><table><thead><tr>' + head +
               '</tr></thead><tbody>' + body + '</tbody></table></html>')
    paths[1].write_text(content, encoding='utf-8')
    paths[2].write_text(json.dumps(dict(summary=summary, residues=rows), indent=2, ensure_ascii=False), encoding='utf-8')
    return [str(path.resolve()) for path in paths]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--input', type=Path, help='FASTA/plain sequence, PDB, mmCIF, or gzip file')
    source.add_argument('--sequence', help='Pasted sequence or FASTA text')
    source.add_argument('--pdb-id', help='RCSB PDB accession')
    source.add_argument('--url', help='Data URL, RCSB entry, UniProt entry, or GPCRdb protein page')
    parser.add_argument('--chain', help='Author or label chain ID; pass an empty string for a blank PDB chain')
    parser.add_argument('--record', help='FASTA record identifier')
    parser.add_argument('--list', action='store_true', help='List records/chains without annotation')
    parser.add_argument('--reference-fasta', type=Path, help='Full sequence for an observed-only structure')
    parser.add_argument('--gpcrdb-entry', help='Verified matching class-A GPCRdb reference, e.g. adrb2_human')
    parser.add_argument('--gpcrdb-cache', type=Path, help='Offline directory containing protein.json and residues.json')
    parser.add_argument('--anchors', default='', help='Evidence-backed full-sequence anchors, e.g. 3.50=131')
    parser.add_argument('--tm5-end', type=int)
    parser.add_argument('--tm6-start', type=int)
    parser.add_argument('--draft', action='store_true', help='Allow explicit draft tables for uncertain inputs')
    parser.add_argument('--out', type=Path, default=Path('gpcr_bw'), help='Output prefix')
    args = parser.parse_args(argv)
    try:
        origin = str(args.input or 'pasted sequence')
        if args.pdb_id:
            if not re.fullmatch(r'(?:[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]{8})', args.pdb_id):
                raise ValueError('Invalid PDB accession.')
            origin = 'https://files.rcsb.org/download/' + args.pdb_id.upper() + '.cif'
            data = fetch(origin)
        elif args.url:
            origin = resolve_url(args.url)
            data = fetch(origin)
        elif args.sequence is not None:
            data = args.sequence.encode()
        else:
            data = args.input.read_bytes()
        records = decode_records(data)
        if args.list:
            print(json.dumps([dict(record=r.name, chain=r.chain, label_chain=r.label_chain,
                                   length=len(r.sequence), observed=len(r.observed), source=r.sequence_source)
                              for r in records], indent=2))
            return 0
        protein = choose(records, args.chain, args.record)
        rows, summary, doubtful = build_table(protein, args)
        summary.update(input_source=origin, input_sha256=hashlib.sha256(data).hexdigest(),
                       created_at_utc=datetime.now(timezone.utc).isoformat())
        if doubtful and not args.draft:
            print(json.dumps(summary, indent=2, ensure_ascii=False))
            print('No final table written. Resolve the warnings, or use --draft for an explicitly provisional table.', file=sys.stderr)
            return 2
        paths = export(rows, summary, args.out, doubtful or args.draft)
        print(json.dumps(dict(summary=summary, files=paths), indent=2, ensure_ascii=False))
        return 0
    except (ValueError, KeyError, OSError, IndexError) as exc:
        print('Error: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())

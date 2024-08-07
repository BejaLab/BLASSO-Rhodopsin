from Bio import SeqIO
from bioformats import read_hhr
from collections import defaultdict

aln_file = str(snakemake.input['aln'])
fas_file = str(snakemake.input['fas'])
a2m_files = snakemake.input['a2m']
hhr_files = snakemake.input['hhr']

output_file = str(snakemake.output)

ref_template = snakemake.params['ref']

aln = SeqIO.to_dict(SeqIO.parse(aln_file, 'clustal'))
gaps = set()
template_pos = defaultdict(dict)
for template, record in aln.items():
    q_pos = 0
    t_pos = 0
    for q, t in zip(record.seq, aln[ref_template].seq):
        q_gap = q == '-'
        t_gap = t == '-'
        q_pos += not q_gap
        t_pos += not t_gap
        if not q_gap and not t_gap:
            template_pos[template][q_pos] = t_pos
        if q_gap and not t_gap:
            gaps.add(t_pos)

best_scores = {}
best_matches = {}

template_seqs = {}
for a2m_file, hhr_file in zip(a2m_files, hhr_files):
    template = None
    a2m_records = {}
    for record in SeqIO.parse(a2m_file, 'fasta'):
        if template is None:
            template = record.id.rsplit('_consensus', maxsplit = 1)[0]
            template_seqs[template] = record.seq
        else:
            rec_num, record.id, *rest = record.description.split(' ') # revert back to the original ename
            record.description = ' '.join([ record.id ] + rest)
            a2m_records[rec_num] = record

    with open(hhr_file) as hhr:
        for record in read_hhr(hhr):
            rec_num = record['query']['name']
            name = a2m_records[rec_num].id
            score = record['Score']
            if name not in best_scores or score > best_scores[name]:
                best_scores[name] = score
                best_matches[name] = a2m_records[rec_num], template

def write_fasta(name, seq, file):
    file.write(f">{name}\n{seq}\n")

with open(output_file, 'w') as out:
    seq = ''
    pos = 0
    for res in aln[ref_template].seq:
        if res != '-':
            pos += 1
            if pos not in gaps:
                seq += res
    write_fasta(ref_template + '_consensus', seq, out)
    for record in SeqIO.parse(fas_file, 'fasta'):
        name = record.id
        a2m_record, template = best_matches[name]
        seq = ''
        t_pos = 0
        for q, t in zip(a2m_record.seq, template_seqs[template]):
            t_gap = t == '.'
            t_pos += not t_gap
            if q != '.':
                ref_pos = template_pos[template][t_pos] if t_pos in template_pos[template] else None
                no_match = ref_pos is None or ref_pos in gaps
                if not (q == '-' and no_match):
                    seq += q if not no_match else q.lower()
        write_fasta(record.description, seq, out)

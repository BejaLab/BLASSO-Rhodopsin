from bioformats import read_hhr
from Bio import SeqIO

fasta_file = str(snakemake.input['fasta'])
hhm_file = str(snakemake.input['hhm'])
hhr_file = str(snakemake.input['hhr'])
output_file = str(snakemake.output)

fasta_records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

def get_cons_from_hhm(file):
    cons_name = None
    cons_seq = ''
    seq_block = False
    for line in file:
        line = line.strip()
        if line == 'SEQ':
            seq_block = True
        elif line == '#':
            seq_block = False
        elif seq_block:
            if line.startswith('>'):
                cons_name = line[1:] if line.endswith('_consensus') else None
            elif cons_name is not None:
                cons_seq += line
    return cons_name, cons_seq

def write_fasta(name, seq, file):
    file.write(f">{name}\n{seq}\n")

with open(hhm_file) as hhm:
    cons_name, cons_seq = get_cons_from_hhm(hhm)

with open(hhr_file) as hhr, open(output_file, 'w') as out:
    write_fasta(cons_name, cons_seq, out)
    for record in read_hhr(hhr):
        if record['Aligned_cols'] > 1:
            query = record['query']
            template = record['template']
            fasta_record = fasta_records[query['name']]
            q_start, q_end = query['coords']
            t_start, t_end = template['coords']

            seq = fasta_record.seq[:q_start - 1].lower()
            seq += '-' * (t_start - 1)
            for q, t in zip(query['alignment'], template['alignment']):
                seq += q if t != '-' else q.lower()
            seq += '-' * (template['len'] - t_end)
            seq += fasta_record.seq[q_end:].lower()

            write_fasta(fasta_record.description, seq, out)

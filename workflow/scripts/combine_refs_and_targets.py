from Bio import SeqIO
from Bio.Seq import Seq

ref_fas_file = str(snakemake.input['ref'])
a3m_file = str(snakemake.input['a3m'])
out_fas_file = str(snakemake.output['fas'])
out_pos_file = str(snakemake.output['pos'])

with open(out_fas_file, 'w') as out_fas, open(out_pos_file, 'w') as out_pos:
    fas_records = SeqIO.parse(ref_fas_file, 'fasta')
    template = next(fas_records)
    for record in fas_records:
        SeqIO.write(record, out_fas, 'fasta')
    a3m_records = SeqIO.parse(a3m_file, 'fasta')
    a3m_template = next(a3m_records)
    assert template.seq.replace('-', '') == a3m_template.seq, "Templates have different sequences"
    for record in a3m_records:
        q_pos = 0
        t_pos = 0
        a3m_mapping = {}
        seq = ''
        for q in record.seq:
            q_gap = q == '-'
            q_pos += not q_gap
            if q == q.upper():
                t_pos += 1
                seq += q
                a3m_mapping[t_pos] = q_pos if not q_gap else -1
        fas_mapping = {}
        t_pos = 0
        gapped_seq = ''
        for i, t in enumerate(template):
            if t == '-':
                gapped_seq += '-'
            else:
                t_pos += 1
                gapped_seq += seq[t_pos - 1]
                q_pos = a3m_mapping[t_pos]
                fas_mapping[q_pos] = i + 1
        
        out_pos.write(f'>{record.description}\n')
        q_pos = 0
        for q in record.seq:
            if q != '-':
                q_pos += 1
                fas_pos = fas_mapping[q_pos] if q_pos in fas_mapping else '-'
                out_pos.write(f"{q}, {q_pos}, {fas_pos}\n")

        record.seq = Seq(gapped_seq)
        SeqIO.write(record, out_fas, 'fasta')

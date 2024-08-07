from Bio import SeqIO
from collections import defaultdict, Counter

aln1_file = snakemake.input['aln1']
aln2_file = snakemake.input['aln2']
pos_input_file = snakemake.input['pos']

pos_output_file = snakemake.output['pos']
counts_file = snakemake.output['counts']

with open(pos_input_file) as fp:
    pocket = [ int(line) for line in fp ]

def aln2seq(seq, aln_positions):
    seq_positions = {}
    seq_pos = 0
    for aln_i, res in enumerate(seq):
        is_not_gap = res != '-'
        seq_pos += is_not_gap
        if aln_i + 1 in aln_positions and is_not_gap:
            seq_positions[seq_pos] = aln_i + 1
    return seq_positions

def seq2aln(seq, seq_positions):
    aln_positions = {}
    seq_pos = 0
    for aln_i, res in enumerate(seq):
        is_not_gap = res != '-'
        seq_pos += is_not_gap
        if seq_pos in seq_positions and is_not_gap:
            aln_positions[aln_i + 1] = seq_pos
    return aln_positions

seq_pocket = {} 
for record in SeqIO.parse(aln1_file, 'fasta'):
    seq_pocket[record.id] = aln2seq(record.seq, pocket)

aln_pocket = defaultdict(list)
for record in SeqIO.parse(aln2_file, 'fasta'):
    if record.id in seq_pocket:
        record_pocket = seq_pocket[record.id]
        for aln2_pos, seq_pos in seq2aln(record.seq, record_pocket).items():
            aln1_pos = record_pocket[seq_pos]
            aln_pocket[aln1_pos].append(aln2_pos)

with open(pos_output_file, 'w') as pos_fp, open(counts_file, 'w') as counts_fp:
    for aln1_pos, aln2_positions in aln_pocket.items():
        counts = Counter(aln2_positions)
        for aln2_pos, count in counts.most_common():
            counts_fp.write(f"{aln1_pos}\t{aln2_pos}\t{count}\n")
        most_common_pos = counts.most_common()[0][0]
        pos_fp.write(f"{most_common_pos}\n")

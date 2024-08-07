from collections import defaultdict
from Bio import SeqIO

tsv_file = snakemake.input['tsv']
fasta_file = snakemake.input['fasta']
dat_file = snakemake.output['dat']
idx_file = snakemake.output['idx']

parse = SeqIO.parse(fasta_file, "fasta")
records = { record.id: record for record in parse }

reprs = defaultdict(list)
with open(tsv_file) as file:
    for line in file:
        cluster_repr, seq_name = line.split()
        reprs[cluster_repr].append(seq_name)

with open(dat_file, 'wb') as dat, open(idx_file, 'w') as idx:
    current = 0
    item = 0
    for cluster_repr, seq_names in reprs.items():
        start = dat.tell()
        for seq_name in seq_names:
            record = records[seq_name]
            rec_str = f'>{record.id}\n{record.seq}\n'
            dat.write(rec_str.encode('utf-8'))
        dat.write(b'\x00')
        start = current
        current = dat.tell()
        length = current - start
        idx.write(f'Cluster{item}\t{start}\t{length}\n')
        item += 1

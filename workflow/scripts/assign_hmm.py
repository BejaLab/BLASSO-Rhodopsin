from bioformats import read_hmmsearch
from Bio import SeqIO
hmmsearch_files = snakemake.input['hmmsearch']
hmmalign_files = snakemake.input['hmmalign']
output_files = snakemake.output

records = {}
for i, hmmsearch_file in enumerate(hmmsearch_files):
    with open(hmmsearch_file) as file:
        for record in read_hmmsearch(file):
            seq_name = record["seq_name"]
            record['profile'] = i
            if seq_name not in records or record['full']['score'] > records[seq_name]['full']['score']:
                records[seq_name] = record

for i, hmmalign_file in enumerate(hmmalign_files):
    with open(output_files[i], 'w') as file:
        for seq_record in SeqIO.parse(hmmalign_file, 'fasta'):
            if records[seq_record.id]['profile'] == i:
                SeqIO.write(seq_record, file, 'fasta')

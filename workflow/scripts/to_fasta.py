from Bio import SeqIO
import warnings

pdb_files = snakemake.input
fas_file = str(snakemake.output['fas'])
txt_file = str(snakemake.output['txt'])
names = snakemake.params['names']

with open(fas_file, 'w') as fas, open(txt_file, 'w') as txt:
    for name, pdb_file in zip(names, pdb_files):
        pdb = SeqIO.parse(pdb_file, "pdb-atom")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            record = next(pdb)
            record.id = name
            record.description = ''
        SeqIO.write(record, fas, "fasta")
        txt.write(">%s  _P_ %s\n" % (record.id, pdb_file))

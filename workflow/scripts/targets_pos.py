
mapout_file = snakemake.input['mapout']
positions_file = snakemake.input['positions']
output_file = str(snakemake.output)

with open(positions_file) as fp:
    pos = { int(line): i for i, line in enumerate(fp) }

with open(output_file, 'w') as output, open(mapout_file) as mapout:
    output.write("name\tsequence_pos\tfasta_pos\tpocket_pos\n")
    for line in mapout:
        if line.startswith('>'):
            name, *description = line[1:].split()
        elif not line.startswith('#'):
            res, orig_pos, aln_pos = line.split(',')
            if aln_pos.strip() != '-':
                res = res.strip()
                orig_pos = int(orig_pos)
                aln_pos = int(aln_pos)
                if aln_pos in pos:
                    output.write(f"{name}\t{orig_pos}\t{aln_pos}\t{pos[aln_pos]}\n")

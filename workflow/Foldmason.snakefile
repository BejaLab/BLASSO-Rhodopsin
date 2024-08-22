rule transfer_positions1:
    input:
        pos = "training/original/position_retinal.dat",
        aln1 = "training/original/sequences.fas",
        aln2 = "training/original/sequences.mafft"
    output:
        pos = "mafft.dat",
        counts = "mafft.counts"
    conda:
        "envs/python.yaml"
    script:
        "scripts/transfer_positions.py"

rule transfer_positions2:
    input:
        pos = "training/original/position_retinal.dat",
        aln1 = "training/original/sequences.fas",
        aln2 = "foldmason_renamed.fas"
    output:
        pos = "foldmason.dat",
        counts = "foldmason.counts"
    conda:
        "envs/python.yaml"
    script:
        "scripts/transfer_positions.py"


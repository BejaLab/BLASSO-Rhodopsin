
ref_template = '7TM'

rule uniref_clust:
    input:
        "resources/UniRef100_domains.faa"
    output:
        "analysis/uniref/uniref_cluster.tsv"
    params:
        db = lambda w, output: output[0].replace("_cluster.tsv", ""),
        min_seq_id = 0.6,
        cov = 0.95,
        cov_mode = 1,
        clu_mode = 2
    conda:
        "envs/mmseqs2.yaml"
    shadow:
        "minimal"
    threads:
        20
    shell:
        "mmseqs easy-cluster {input} {params.db} tmp --threads {threads} --min-seq-id {params.min_seq_id} -c {params.cov} --cov-mode {params.cov_mode} --cluster-mode {params.clu_mode}"

rule uniref_ffindex:
    input:
        tsv = "analysis/uniref/uniref_cluster.tsv",
        fasta = "resources/UniRef100_domains.faa"
    output:
        dat = "analysis/uniref/uniref_faa.ffdata",
        idx = "analysis/uniref/uniref_faa.ffindex"
    conda:
        "envs/python.yaml"
    script:
        "scripts/cluster_ffindex.py"

rule uniref_mafft:
    input:
        dat = "{prefix}_faa.ffdata",
        idx = "{prefix}_faa.ffindex"
    output:
        dat = "{prefix}_msa.ffdata",
        idx = "{prefix}_msa.ffindex"
    log:
        "{prefix}_msa.log"
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        "mpirun --oversubscribe -np {threads} ffindex_apply_mpi -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- mafft --auto - &> {log}"

rule uniref_a3m:
    input:
        dat = "{prefix}_msa.ffdata",
        idx = "{prefix}_msa.ffindex"
    output:
        dat = "{prefix}_msa_a3m.ffdata",
        idx = "{prefix}_msa_a3m.ffindex"
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        mpirun --oversubscribe -np {threads} ffindex_apply_mpi -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- reformat.pl -v 0 -M first fas a3m /dev/stdin /dev/stdout
        """

rule uniref_cstranslate_a3m:
    input:
        dat = "{prefix}_msa_a3m.ffdata",
        idx = "{prefix}_msa_a3m.ffindex"
    output:
        dat = "{prefix}_msa_cs219.ffdata",
        idx = "{prefix}_msa_cs219.ffindex"
    log:
        "{prefix}_msa_cs219.log"
    params:
        DB_in  = "{prefix}_msa_a3m",
        DB_out = "{prefix}_msa_cs219"
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        export DATA_DIR=$CONDA_PREFIX/data
        mpirun --oversubscribe -np {threads} cstranslate_mpi -A $DATA_DIR/cs219.lib -D $DATA_DIR/context_data.lib -i {params.DB_in} -o {params.DB_out} -x 0.3 -c 4 -I a3m -b &> {log}
        """
profiles ,= glob_wildcards("profiles/{profile}.hhm")

rule copy_pdb:
    input:
        "profiles/{profile}_colabfold.pdb"
    output:
        "analysis/profiles/{profile}.pdb"
    shell:
        "cp {input} {output}"

rule pdb_to_fasta:
    input:
        expand("analysis/profiles/{profile}.pdb", profile = profiles)
    output:
        fas = "analysis/profiles/sequences.fasta",
        txt = "analysis/profiles/template.txt"
    params:
        names = profiles
    conda:
        "envs/python.yaml"
    script:
        "scripts/to_fasta.py"

rule t_coffee:
    input:
        fas = "analysis/profiles/sequences.fasta",
        txt = "analysis/profiles/template.txt"
    output:
        "analysis/profiles/alignment.aln"
    params:
        method = "sap_pair,mustang_pair"
    log:
        "analysis/profiles/alignment.log"
    conda:
        "envs/t_coffee.yaml"
    shell:
        "t_coffee {input.fas} -outfile {output} -method {params.method} -template_file {input.txt} -n_core 1 &> {log}"

rule numbers_as_names:
    input:
        "{prefix}/sequences.fasta"
    output:
        "{prefix}/sequences_num.fasta"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit fx2tab {input} | nl -w1 | sed 's/\\t/ /' | seqkit tab2fx -o {output}"

rule ffindex_from_fasta:
    input:
        "{prefix}/sequences_num.fasta"
    output:
        dat = "{prefix}/sequences.ffdata",
        idx = "{prefix}/sequences.ffindex"
    conda:
        "envs/hhsuite.yaml"
    shell:
        "ffindex_from_fasta -s {output.dat} {output.idx} {input}"

rule uniref_hhblits:
    input:
        dat = "{prefix}/sequences.ffdata",
        idx = "{prefix}/sequences.ffindex",
        DB = expand("analysis/uniref/uniref_msa_{type}.{ext}", type = [ "a3m", "cs219" ], ext = [ "ffdata", "ffindex" ])
    output:
        dat = "{prefix}/sequences_a3m.ffdata",
        idx = "{prefix}/sequences_a3m.ffindex"
    log:
        "{prefix}/sequences_a3m.log"
    params:
        db = "analysis/uniref/uniref_msa",
        e = 1e-5,
        n = 2
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        mpirun --oversubscribe -np {threads} ffindex_apply_mpi -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- \
            hhblits -cpu 1 -v 0 -b 1 -z 1 -d {params.db} -i stdin -oa3m stdout -o /dev/null -e {params.e} -n {params.n} &> {log}
        """

rule addss:
    input:
        dat = "{prefix}/sequences_a3m.ffdata",
        idx = "{prefix}/sequences_a3m.ffindex"
    output:
        dat = "{prefix}/sequences_ss_a3m.ffdata",
        idx = "{prefix}/sequences_ss_a3m.ffindex"
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        mpirun --oversubscribe -np {threads} ffindex_apply_mpi -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- addss.pl -a3m -v 0 stdin stdout
        """
rule hhalign:
    input:
        dat = "{prefix}/sequences_a3m.ffdata",
        idx = "{prefix}/sequences_a3m.ffindex",
        hhm = "profiles/{profile}.hhm"
    output:
        dat = "{prefix}/hhalign/{profile}_hhr.ffdata",
        idx = "{prefix}/hhalign/{profile}_hhr.ffindex"
    log:
        "{prefix}/hhalign/{profile}_hhr.log"
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        mpirun --oversubscribe -np {threads} ffindex_apply_mpi -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- \
            hhalign -v 0 -glob -i stdin -t {input.hhm} -o stdout &> {log}
        """

rule hhr_to_a3m:
    input:
        hhm = "profiles/{profile}.hhm",
        hhr = "{prefix}/hhalign/{profile}_hhr.ffdata",
        fasta = "{prefix}/sequences_num.fasta"
    output:
        "{prefix}/hhalign/{profile}.a3m"
    conda:
        "envs/python.yaml"
    script:
        "scripts/hhr_to_a3m.py"

rule a3m_to_a2m:
    input:
        "{prefix}/hhalign/{profile}.a3m"
    output:
        "{prefix}/hhalign/{profile}.a2m"
    conda:
        "envs/hhsuite.yaml"
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        reformat.pl a3m a2m {input} {output}
        """

rule aggregate_alignments:
    input:
        aln = "analysis/profiles/alignment.aln",
        fas = "{prefix}/sequences.fasta",
        a2m = expand("{{prefix}}/hhalign/{profile}.a2m", profile = profiles),
        hhr = expand("{{prefix}}/hhalign/{profile}_hhr.ffdata", profile = profiles)
    output:
        "{prefix}/hhalign/profiles.a3m"
    params:
        ref = ref_template
    conda:
        "envs/python.yaml"
    script:
        "scripts/aggregate_alignments.py"

rule a3m_to_fas:
    input:
        "{prefix}.a3m"
    output:
        "{prefix}_a3m.fas"
    conda:
        "envs/hhsuite.yaml"
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        reformat.pl a3m fas {input} {output}
        """

rule profile_ref_aln:
    input:
        "training/{set}-profiles/hhalign/profiles_a3m.fas"
    output:
        "training/{set}-profiles/sequences.fas"
    params:
        ref = ref_template
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit grep -vp {params.ref}_consensus {input} > {output}"

rule transfer_positions:
    input:
        pos = "training/{set}/position_retinal.dat",
        aln1 = "training/{set}/sequences.fas",
        aln2 = "training/{set}-profiles/sequences.fas"
    output:
        pos = "training/{set}-profiles/position_retinal.dat",
        counts = "training/{set}-profiles/position_retinal.counts"
    conda:
        "envs/python.yaml"
    script:
        "scripts/transfer_positions.py"

rule combine_refs_and_targets:
    input:
        ref = "training/{set}-profiles/hhalign/profiles_a3m.fas",
        a3m = "analysis/targets/{target}/hhalign/profiles.a3m"
    output:
        fas = "analysis/predict/{set}-profiles/{target}.fas",
        pos = "analysis/predict/{set}-profiles/{target}.txt"
    conda:
        "envs/python.yaml"
    script:
        "scripts/combine_refs_and_targets.py"

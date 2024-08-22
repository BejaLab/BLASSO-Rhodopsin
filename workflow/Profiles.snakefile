
ref_template = '7TM'
profiles ,= glob_wildcards("resources/profiles/{profile}.hhm")

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
        "ffindex_apply_py -j {threads} -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- mafft --auto - &> {log}"

rule uniref_a3m:
    input:
        dat = "{prefix}_msa.ffdata",
        idx = "{prefix}_msa.ffindex"
    output:
        dat = "{prefix}_msa_a3m.ffdata",
        idx = "{prefix}_msa_a3m.ffindex"
    log:
        "{prefix}_msa_a3m.log"
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        ffindex_apply_py -j {threads} -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- reformat.pl -v 0 -M first fas a3m /dev/stdin /dev/stdout &> {log}
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
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        export DATA_DIR=$CONDA_PREFIX/data
        ffindex_apply_py -j {threads}  -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- \
            cstranslate -A $DATA_DIR/cs219.lib -D $DATA_DIR/context_data.lib -i stdin -o stdout -x 0.3 -c 4 -I a3m -b &> {log}
        """

rule ffindex_from_fasta:
    input:
        "{prefix}/sequences.fasta"
    output:
        dat = "{prefix}/sequences.ffdata",
        idx = "{prefix}/sequences.ffindex"
    conda:
        "envs/hhsuite.yaml"
    shell:
        "ffindex_from_fasta_py {output.dat} {output.idx} {input}"

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
        qid = 30,
        e = 1e-10,
        n = 2
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        ffindex_apply_py -j {threads} -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- \
            hhblits -cpu 1 -v 0 -b 1 -z 1 -d {params.db} -i stdin -oa3m stdout -o /dev/null -e {params.e} -n {params.n} -qid {params.qid} &> {log}
        """

rule addss:
    input:
        dat = "{prefix}/sequences_a3m.ffdata",
        idx = "{prefix}/sequences_a3m.ffindex"
    output:
        dat = "{prefix}/sequences_ss_a3m.ffdata",
        idx = "{prefix}/sequences_ss_a3m.ffindex"
    log:
        "{prefix}/sequences_ss_a3m.log"
    conda:
        "envs/hhsuite.yaml"
    threads:
        workflow.cores
    shell:
        """
        export PATH=$CONDA_PREFIX/scripts:$PATH
        ffindex_apply_py -j {threads} -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- addss.pl -a3m -v 0 stdin stdout &> {log}
        """

rule hhalign:
    input:
        dat = "{prefix}/sequences_a3m.ffdata",
        idx = "{prefix}/sequences_a3m.ffindex",
        hhm = "resources/profiles/{profile}.hhm"
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
        "ffindex_apply_py -j {threads} -d {output.dat} -i {output.idx} {input.dat} {input.idx} -- hhalign -v 0 -glob -i stdin -t {input.hhm} -o stdout &> {log}"

rule hhr_to_a3m:
    input:
        hhm = "resources/profiles/{profile}.hhm",
        hhr = "{prefix}/hhalign/{profile}_hhr.ffdata",
        fasta = "{prefix}/sequences.fasta"
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
        aln = "resources/profiles.aln",
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
        "envs/seqkit.yaml"
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

rule targets_pos_profiles:
    input:
        mapout = "analysis/predict/{set}-profiles/{target}.txt",
        positions = "training/{set}-profiles/position_retinal.dat"
    output:
        "output/{set}-profiles/{target}.pos"
    conda:
        "envs/python.yaml"
    script:
        "scripts/targets_pos.py"

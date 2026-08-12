# prepare samplesheet
# -----------------------------------------------------
rule samplesheet:
    input:
        path=config["samplesheet"],
    output:
        path="results/samplesheet/samplesheet.tsv",
    log:
        path="results/samplesheet/samplesheet.log",
    conda:
        "../envs/basic.yml"
    script:
        "../scripts/prepare_samplesheet.py"


# fetch protein database from NCBI
# -----------------------------------------------------
rule database:
    input:
        fasta=config["database"],
    output:
        fasta="results/database/database.fasta",
    log:
        path="results/database/database.log",
    conda:
        "../envs/basic.yml"
    script:
        "../scripts/prepare_database.py"


# generate decoys
# -----------------------------------------------------
rule decoypyrat:
    input:
        path=rules.database.output.fasta,
    output:
        path="results/decoypyrat/decoy_database.fasta",
    log:
        path="results/decoypyrat/decoypyrat.log",
    conda:
        "../envs/decoypyrat.yml"
    params:
        cleavage_sites=config["decoypyrat"]["cleavage_sites"],
        decoy_prefix=config["decoypyrat"]["decoy_prefix"],
    shell:
        "if ! grep -q '>rev_' {input.path};" "then decoypyrat {input.path} \
        -c {params.cleavage_sites} \
        -d {params.decoy_prefix} \
        -o {output.path} \
        -k > {log.path}; fi;" "cat {input.path} >> {output.path}"


# prepare workflow
# -----------------------------------------------------
rule workflow:
    input:
        samplesheet=rules.samplesheet.output.path,
        database=rules.decoypyrat.output.path,
    output:
        path="results/workflow/workflow.txt",
    log:
        path="results/workflow/workflow.log",
    conda:
        "../envs/basic.yml"
    params:
        workflow=config["workflow"],
        workflow_dir=wfpath("../resources"),
    script:
        "../scripts/prepare_workflow.py"

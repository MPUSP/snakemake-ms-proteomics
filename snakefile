# ----------------------------------------------------- #
# snakemake pipeline for automatic analysis/QC of mass  #
# spectrometry data                                     #
#                                                       #
# Author: Michael Jahn                                  #
# Date: 2022-11-24                                      #
# License: GPL v3 (for all 3rd party tools              #
# separate licenses apply)                              #
# ----------------------------------------------------- #

# import basic packages
from os import path


# check command line arguments
# -----------------------------------------------------
print("\n +++ SNAKEMAKE GLOBAL PARAMETERS +++ \n")
for i in config.keys():
    print(f"  - {i}: {config.get(i)}")
print("\n")


# function to construct target paths
def out(file):
    return os.path.join(config["output"], file)


# target rule
# -----------------------------------------------------
rule all:
    input:
        out("email/log.txt"),
        out("clean_up/log.txt"),
        out("versions/log_env.txt"),
        out("module_logs/log.txt"),


# module to fetch protein database from NCBI
# -----------------------------------------------------
rule database:
    params:
        term=config["database"],
    output:
        path=directory(out("database")),
        database=out("database/database.fasta"),
    conda:
        "envs/database.yml"
    log:
        path=out("database/log.txt"),
    script:
        "scripts/prepare_database.py"


# module to generate decoys
# -----------------------------------------------------
rule decoypyrat:
    input:
        path=rules.database.output.database,
    output:
        path=out("decoypyrat/decoy_database.fasta"),
    conda:
        "envs/decoypyrat.yml"
    params:
        cleavage_sites=config["decoypyrat"]["cleavage_sites"],
        decoy_prefix=config["decoypyrat"]["decoy_prefix"],
    log:
        path=out("decoypyrat/log.txt"),
    shell:
        "if ! grep -q '>rev_' {input.path};"
        "then decoypyrat {input.path} \
        -c {params.cleavage_sites} \
        -d {params.decoy_prefix} \
        -o {output.path} \
        -k > {log.path}; fi;"
        "cat {input.path} >> {output.path}"


# module to prepare samplesheet
# -----------------------------------------------------
rule samplesheet:
    input:
        path=config["samplesheet"],
    output:
        path=out("samplesheet/samplesheet.tsv"),
    conda:
        "envs/samplesheet.yml"
    log:
        path=out("samplesheet/log.txt"),
    script:
        "scripts/prepare_samplesheet.py"


# module to prepare workflow
# -----------------------------------------------------
rule workflow:
    input:
        samplesheet=rules.samplesheet.output.path,
        database=rules.decoypyrat.output.path,
    output:
        path=out("workflow/workflow.txt"),
    conda:
        "envs/workflow.yml"
    params:
        workflow=config["workflow"],
    log:
        path=out("workflow/log.txt"),
    script:
        "scripts/prepare_workflow.py"


# module to run fragpipe
# -----------------------------------------------------
rule fragpipe:
    input:
        fragpipe_bin=config["fragpipe"]["path"],
        samplesheet=rules.samplesheet.output.path,
        workflow=rules.workflow.output.path,
    output:
        path=directory(out("fragpipe")),
        msstats=out("fragpipe/MSstats.csv"),
    conda:
        "envs/fragpipe.yml"
    params:
        dummyParam=0,
    log:
        path=out("fragpipe/log.txt"),
    shell:
        "{input.fragpipe_bin}/fragpipe \
        --headless \
        --workflow {input.workflow} \
        --manifest {input.samplesheet} \
        --workdir {output.path} \
        > {log.path}"


# module to run MSstats
# -----------------------------------------------------
rule msstats:
    input:
        samplesheet=rules.samplesheet.output.path,
        table_msstats=rules.fragpipe.output.msstats,
    output:
        feature_level_data=out("msstats/feature_level_data.csv"),
        protein_level_data=out("msstats/protein_level_data.csv"),
        comparison_result=out("msstats/comparison_result.csv"),
        model_qc=out("msstats/model_qc.csv"),
        uniprot=out("msstats/uniprot.csv"),
    conda:
        "envs/msstats.yml"
    params:
        config_msstats=config["msstats"],
    log:
        path=out("msstats/log.txt"),
    script:
        "scripts/run_msstats.R"


# module to clean up files after pipeline execution
# -----------------------------------------------------
rule clean_up:
    input:
        samplesheet=rules.samplesheet.output.path,
        msstats=rules.fragpipe.output.msstats,
    output:
        log=out("clean_up/log.txt"),
    params:
        pattern="_uncalibrated.mzML",
    shell:
        "echo 'removed the following files:' >> {output.log};"
        "while read -r line;"
        "do filename=`echo ${{line}} | cut -f 1 -d ' '`;"
        "filename=`echo ${{filename//.raw/{params.pattern}}}`;"
        "if test -f ${{filename}}; then rm ${{filename}}; echo ${{filename}} >> {output.log}; fi;"
        "done < {input.samplesheet};"


# module to fetch software versions from conda envs
# -----------------------------------------------------
rule versions:
    input:
        expand(
            "envs/{module}.yml",
            module=[
                "database",
                "decoypyrat",
                "email",
                "fragpipe",
                "msstats",
                "pdf",
                "report",
                "samplesheet",
                "workflow",
            ],
        ),
    output:
        path=out("versions/log_packages.txt"),
    log:
        path=out("versions/log_env.txt"),
    shell:
        "conda env export > {log.path};"
        "cat {input} >> {output.path}"


# module to combine all module log files to single log
# -----------------------------------------------------
rule module_logs:
    input:
        rules.database.log.path,
        rules.decoypyrat.log.path,
        rules.samplesheet.log.path,
        rules.workflow.log.path,
        rules.fragpipe.log.path,
        rules.msstats.log.path,
        rules.clean_up.output.log,
    log:
        path=out("module_logs/log.txt"),
    shell:
        "cat {input} >> {log.path}"


# module to generate full HTML report using R markdown
# -----------------------------------------------------
rule report:
    input:
        feature_level_data=rules.msstats.output.feature_level_data,
        protein_level_data=rules.msstats.output.protein_level_data,
        comparison_result=rules.msstats.output.comparison_result,
        model_qc=rules.msstats.output.model_qc,
        versions=rules.versions.output.path,
    output:
        html=out("report/report.html"),
    conda:
        "envs/report.yml"
    params:
        config_report=config["report"],
    script:
        "notebooks/report.Rmd"


# module to convert HTML to PDF output
# -----------------------------------------------------
rule pdf:
    input:
        html=rules.report.output.html,
    output:
        pdf=out("report/report.pdf"),
    conda:
        "envs/pdf.yml"
    log:
        path=out("report/log.txt"),
    shell:
        "weasyprint -v {input.html} {output.pdf} &> {log.path}"


# module to send out emails using custom mail server
# -----------------------------------------------------
rule email:
    input:
        html=rules.report.output.html,
        pdf=rules.pdf.output.pdf,
        protein=rules.msstats.output.protein_level_data,
        comparison=rules.msstats.output.comparison_result,
    output:
        log=out("email/log.txt"),
    params:
        config_email=config["email"],
        config_database=config["database"],
        config_workflow=config["workflow"],
        config_samplesheet=config["samplesheet"],
        config_out_dir=config["output"],
    script:
        "scripts/send_email.py"

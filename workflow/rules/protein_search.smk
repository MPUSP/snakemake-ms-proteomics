# download and test fragpipe
# -----------------------------------------------------
rule fragpipe_setup:
    output:
        executable=f"results/fragpipe_setup/{config['fragpipe']['executable']}",
    log:
        path="results/fragpipe_setup/setup.log",
    conda:
        "../envs/fragpipe.yml"
    params:
        fragpipe_download=config["fragpipe"]["download"],
    shell:
        "set -euo pipefail;"
        "if test -f {output.executable}; then exit 0; fi;"
        "wget -O results/fragpipe_setup/fragpipe.zip {params.fragpipe_download} > {log.path} 2>&1;"
        "unzip -o -d results/fragpipe_setup/ results/fragpipe_setup/fragpipe.zip > {log.path} 2>&1;"
        "rm -f results/fragpipe_setup/fragpipe.zip;"
        "test -f {output.executable};"
        "{output.executable}"


# run fragpipe
# -----------------------------------------------------
rule fragpipe:
    input:
        samplesheet=rules.samplesheet.output.path,
        workflow=rules.workflow.output.path,
        executable=rules.fragpipe_setup.output.executable,
    output:
        path=directory("results/fragpipe"),
        msstats="results/fragpipe/msstats.csv",
    log:
        path="results/fragpipe/fragpipe_module.log",
    conda:
        "../envs/fragpipe.yml"
    shell:
        "set -euo pipefail;"
        "{input.executable} "
        "--headless "
        "--workflow {input.workflow} "
        "--manifest {input.samplesheet} "
        "--workdir {output.path} "
        "> {log.path};"
        "if test -f {output.path}/dia-quant-output/msstats.csv;"
        "then cp {output.path}/dia-quant-output/msstats.csv {output.msstats}; fi;"


# run MSstats
# -----------------------------------------------------
rule msstats:
    input:
        samplesheet=rules.samplesheet.output.path,
        table_msstats=rules.fragpipe.output.msstats,
    output:
        feature_level_data="results/msstats/feature_level_data.csv",
        protein_level_data="results/msstats/protein_level_data.csv",
        comparison_result="results/msstats/comparison_result.csv",
        model_qc="results/msstats/model_qc.csv",
        uniprot="results/msstats/uniprot.csv",
    log:
        path="results/msstats/msstats.log",
    conda:
        "../envs/msstats.yml"
    params:
        config_msstats=config["msstats"],
    script:
        "../scripts/run_msstats.R"

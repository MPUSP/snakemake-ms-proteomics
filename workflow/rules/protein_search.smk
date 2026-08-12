# run fragpipe
# -----------------------------------------------------
rule fragpipe:
    input:
        samplesheet=rules.samplesheet.output.path,
        workflow=rules.workflow.output.path,
    output:
        path=directory("results/fragpipe"),
        msstats="results/fragpipe/msstats.csv",
    log:
        path="results/fragpipe/fragpipe_module.log",
    conda:
        "../envs/fragpipe.yml"
    params:
        fragpipe_dir=config["fragpipe"]["target_dir"],
        fragpipe_bin=config["fragpipe"]["executable"],
        fragpipe_download=config["fragpipe"]["download"],
    shell:
        "env=`echo $CONDA_PREFIX`;"
        "if ! test -f ${{env}}/{params.fragpipe_dir}/{params.fragpipe_bin};"
        "then wget -P ${{env}}/{params.fragpipe_dir} {params.fragpipe_download};"
        "unzip -d ${{env}}/{params.fragpipe_dir} ${{env}}/{params.fragpipe_dir}/*.zip;"
        "${{env}}/{params.fragpipe_dir}/{params.fragpipe_bin};"
        "fi;"
        "${{env}}/{params.fragpipe_dir}/{params.fragpipe_bin} "
        "--headless "
        "--workflow {input.workflow} "
        "--manifest {input.samplesheet} "
        "--workdir {output.path} "
        "> {log.path};"


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

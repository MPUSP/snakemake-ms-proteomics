# clean up files after pipeline execution
# -----------------------------------------------------
rule clean_up:
    input:
        samplesheet=rules.samplesheet.output.path,
        msstats=rules.fragpipe.output.msstats,
    log:
        path="results/clean_up/clean_up.log",
    conda:
        "../envs/basic.yml"
    params:
        pattern="_uncalibrated.mzML",
    shell:
        "echo 'removed the following files:' >> {log.path};"
        "while read -r line;"
        "do filename=`echo ${{line}} | cut -f 1 -d ' '`;"
        "filename=`echo ${{filename//.raw/{params.pattern}}}`;"
        "if test -f ${{filename}}; then rm ${{filename}}; echo ${{filename}} >> {log.path}; fi;"
        "done < {input.samplesheet};"


# fetch software versions from conda envs
# -----------------------------------------------------
rule versions:
    input:
        expand(
            "../workflow/../envs/{module}.yml",
            module=[
                "basic",
                "database",
                "decoypyrat",
                "email",
                "fragpipe",
                "msstats",
                "report_html",
                "report_pdf",
                "workflow",
            ],
        ),
    output:
        path="results/versions/packages.txt",
    log:
        path="results/versions/versions.log",
    conda:
        "../envs/basic.yml"
    shell:
        "conda env export > {log.path};" "cat {input} >> {output.path}"


# combine all module log files to single log
# -----------------------------------------------------
rule module_logs:
    input:
        rules.database.log.path,
        rules.decoypyrat.log.path,
        rules.samplesheet.log.path,
        rules.workflow.log.path,
        rules.fragpipe.log.path,
        rules.msstats.log.path,
        rules.clean_up.log.path,
    log:
        path="results/module_logs/all.log",
    conda:
        "../envs/basic.yml"
    shell:
        "cat {input} >> {log.path}"


# generate full HTML report using R markdown
# -----------------------------------------------------
rule report_html:
    input:
        feature_level_data=rules.msstats.output.feature_level_data,
        protein_level_data=rules.msstats.output.protein_level_data,
        comparison_result=rules.msstats.output.comparison_result,
        model_qc=rules.msstats.output.model_qc,
        versions=rules.versions.output.path,
    output:
        html="results/report/report.html",
    log:
        path="results/report/report_html.log",
    conda:
        "../envs/report_html.yml"
    params:
        config_report=config["report"],
    script:
        "../notebooks/report.Rmd"


# convert HTML to PDF output
# -----------------------------------------------------
rule report_pdf:
    input:
        html=rules.report_html.output.html,
    output:
        pdf="results/report/report.pdf",
    log:
        path="results/report/report_pdf.log",
    conda:
        "../envs/report_pdf.yml"
    shell:
        "weasyprint -v {input.html} {output.pdf} &> {log.path}"


# send out emails using custom mail server
# -----------------------------------------------------
rule email:
    input:
        html=rules.report_html.output.html,
        pdf=rules.report_pdf.output.pdf,
        protein=rules.msstats.output.protein_level_data,
        comparison=rules.msstats.output.comparison_result,
    output:
        path=directory("results/email"),
    log:
        path="results/email/email.log",
    conda:
        "../envs/basic.yml"
    params:
        config_email=config["email"],
        config_database=config["database"],
        config_workflow=config["workflow"],
        config_samplesheet=config["samplesheet"],
    script:
        "../scripts/send_email.py"

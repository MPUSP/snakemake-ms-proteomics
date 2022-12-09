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
    return(os.path.join(config['output'], file))


# target rule
# -----------------------------------------------------
rule all:
    input:
        out('report/report.html'),
        out('clean_up/log.txt')


# module to fetch protein database from NCBI
# -----------------------------------------------------
rule database:
    params:
        term = config['database']
    output:
        path = directory(out('database')),
        database = out('database/database.fasta'),
    script:
        "source/prepare_database.py"


# module to generate decoys
# -----------------------------------------------------
rule decoypyrat:
    input:
        path = rules.database.output.database
    output:
        path = out('decoypyrat/decoy_database.fasta')
    params:
        cleavage_sites = config['decoypyrat']['cleavage_sites'],
        decoy_prefix = config['decoypyrat']['decoy_prefix']
    shell:
        "if ! grep -q '>rev_' {input.path};"
        "then decoypyrat {input.path} \
        -c {params.cleavage_sites} \
        -d {params.decoy_prefix} \
        -o {output.path} \
        -k; fi;"
        "cat {input.path} >> {output.path}"


# module to prepare workflow
# -----------------------------------------------------
rule workflow:
    input:
        samplesheet = config['samplesheet'],
        database = rules.decoypyrat.output.path
    output:
        path = out('workflow/workflow.txt')
    params:
        workflow = config['workflow']
    script:
        "source/prepare_workflow.py"


# module to run fragpipe
# -----------------------------------------------------
rule fragpipe:
    input:
        fragpipe_bin = config['fragpipe']['path'],
        samplesheet = config['samplesheet'],
        workflow = rules.workflow.output.path
    output:
        path = directory(out('fragpipe')),
        msstats = out('fragpipe/MSstats.csv')
    params:
        dummyParam = 0
    shell:
        "{input.fragpipe_bin}/fragpipe \
        --headless \
        --workflow {input.workflow} \
        --manifest {input.samplesheet} \
        --workdir {output.path}"


# module to run MSstats
# -----------------------------------------------------
rule msstats:
    input:
        samplesheet = config['samplesheet'],
        table_msstats = rules.fragpipe.output.msstats
    output:
        feature_level_data = out('msstats/feature_level_data.csv'),
        protein_level_data = out('msstats/protein_level_data.csv'),
        comparison_result = out('msstats/comparison_result.csv'),
        model_qc = out('msstats/model_qc.csv'),
        uniprot = out('msstats/uniprot.csv')
    params:
        config_msstats = config['msstats']
    script:
        "source/run_msstats.R"


# module to clean up files after pipeline execution
# -----------------------------------------------------
rule clean_up:
    input:
        samplesheet = config['samplesheet'],
        msstats = rules.fragpipe.output.msstats
    output:
        log = out('clean_up/log.txt')
    params:
        pattern = '_uncalibrated.mzML'
    shell:
        "echo 'removed the following files:' >> {output.log};"
        "while read -r line;"
        "do filename=`echo ${{line}} | cut -f 1 -d ' '`;"
        "filename=`echo ${{filename//.raw/{params.pattern}}}`;"
        "if test -f ${{filename}}; then rm ${{filename}}; echo ${{filename}} >> {output.log}; fi;"
        "done < {input.samplesheet};"


# module to generate full HTML report using R markdown
# -----------------------------------------------------
rule report:
    input:
        feature_level_data = rules.msstats.output.feature_level_data,
        protein_level_data = rules.msstats.output.protein_level_data,
        comparison_result = rules.msstats.output.comparison_result,
        model_qc = rules.msstats.output.model_qc
    output:
        log = out('report/report.html')
    params:
        config_report = config['report']
    script:
        "source/report.Rmd"

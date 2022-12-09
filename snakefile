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
import os


# check command line arguments
# -----------------------------------------------------
print("\n +++ SNAKEMAKE GLOBAL PARAMETERS +++ \n")
for i in config.keys():
    print(f"  - {i}: {config.get(i)}")
print("\n")


# target rule
# -----------------------------------------------------
rule all:
    input:
        os.path.join(config['output'], 'decoypyrat/decoy_database.fasta')
        #os.path.join(config['output'], 'report/report.html'),
        #os.path.join(config['output'], 'clean_up/log.txt')


# module to fetch protein database from NCBI
# -----------------------------------------------------
rule database:
    params:
        term = config['database']
    output:
        path = directory(os.path.join(config['output'], 'database')),
        database = os.path.join(config['output'], 'database/database.fasta'),
    script:
        "source/prepare_database.py"


# module to generate decoys
# -----------------------------------------------------
rule decoypyrat:
    input:
        path = os.path.join(config['output'], 'database/database.fasta')
    output:
        path = os.path.join(config['output'], 'decoypyrat/decoy_database.fasta')
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
        database = os.path.join(config['output'], 'decoypyrat/decoy_database.fasta')
    output:
        path = os.path.join(config['output'], 'workflow/workflow.txt')
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
        workflow = os.path.join(config['output'], 'workflow/workflow.txt')
    output:
        path = directory(os.path.join(config['output'], 'fragpipe')),
        msstats = os.path.join(config['output'], 'fragpipe/MSstats.csv')
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
        table_msstats = os.path.join(config['output'], 'fragpipe/MSstats.csv')
    output:
        feature_level_data = os.path.join(config['output'], 'msstats/feature_level_data.csv'),
        protein_level_data = os.path.join(config['output'], 'msstats/protein_level_data.csv'),
        comparison_result = os.path.join(config['output'], 'msstats/comparison_result.csv'),
        model_qc = os.path.join(config['output'], 'msstats/model_qc.csv'),
        uniprot = os.path.join(config['output'], 'msstats/uniprot.csv')
    params:
        config_msstats = config['msstats']
    script:
        "source/run_msstats.R"


# module to clean up files after pipeline execution
# -----------------------------------------------------
rule clean_up:
    input:
        samplesheet = config['samplesheet'],
        msstats = os.path.join(config['output'], 'fragpipe/MSstats.csv')
    output:
        log = os.path.join(config['output'], 'clean_up/log.txt')
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
        feature_level_data = os.path.join(config['output'], 'msstats/feature_level_data.csv'),
        protein_level_data = os.path.join(config['output'], 'msstats/protein_level_data.csv'),
        comparison_result = os.path.join(config['output'], 'msstats/comparison_result.csv'),
        model_qc = os.path.join(config['output'], 'msstats/model_qc.csv'),
    output:
        log = os.path.join(config['output'], 'report/report.html')
    params:
        config_report = config['report']
    script:
        "source/report.Rmd"

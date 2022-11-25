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
        os.path.join(config['output'], 'fragpipe/MSstats.csv')


# module to generate decoys
# -----------------------------------------------------
rule decoypyrat:
    input:
        path = config['database']
    output:
        path = os.path.join(config['output'], 'decoypyrat/decoy_database.fasta')
    shell:
        "decoypyrat {input.path} \
        -d 'rev' \
        -o {output.path} \
        -k;"
        "cat {input.path} >> {output.path}"


# module to prepare workflow
# -----------------------------------------------------
rule workflow:
    input:
        path = config['workflow'],
        database = os.path.join(config['output'], 'decoypyrat/decoy_database.fasta')
    output:
        path = os.path.join(config['output'], 'workflow/workflow.txt')
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

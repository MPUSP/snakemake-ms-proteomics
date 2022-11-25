#!/usr/bin/env Rscript

from os import path

wf_path_in = snakemake.input[0]
wf_path_out = snakemake.output[0]
db_path = snakemake.input[1]
db_var = "database.db-path"

# import workflow and add database term
with open(wf_path_in, "r") as wf_file:
    wf = wf_file.read()
    wf = wf + f"\n{db_var}={path.abspath(db_path)}"

# export workflow
with open(wf_path_out, "w") as wf_out:
    wf_out.write(wf)

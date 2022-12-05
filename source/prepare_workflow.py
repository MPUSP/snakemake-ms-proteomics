#!/usr/bin/env Rscript

# FRAGPIPE WORKFLOWS based on sample type
# -----------------------------------------------------------------------------
#
# if all samples are DDA (data dependent acqisition):
#  --> default workflow is LFQ-MBR
#
# if some samples are one of the following:
# DIA = wide window DIA
# GPF-DIA = gas phase fractionation DIA
# DIA-Quant = only quantification (no ident.?)
# DIA-Lib = only for spectral library generation
#  --> default workflow is DIA_SpecLib_Quant
#
# if user supplies workflow:
#  --> no default workflow is used
#
# Notes:
# Runs with DDA, DIA and GPF-DIA will be used from ident. to quant.
# Runs with DIA-Quant will only be used in quantification

from os import path

wf_path_in = snakemake.params["workflow"]
wf_path_out = snakemake.output["path"]
db_path = snakemake.input["database"]
sample_path = snakemake.input["samplesheet"]
db_var = "database.db-path"

# determine if workflow file was supplied
if wf_path_in == "from_samplesheet":
    # if not, determine wf from sample type
    data_type = []
    with open(sample_path, "r") as sample_file:
        for sf in sample_file.readlines():
            data_type += [sf.split("\t")[3]]
    if all([i == "DDA" for i in data_type]):
        local_path = path.abspath("workflows/LFQ-MBR.workflow")
        print(
            f"WORKFLOW: detected DDA samples, choosing default workflow: {local_path}"
        )
    elif any(["DIA" in i for i in data_type]):
        raise ValueError(
            "WORKFLOW: detected DIA samples, a default workflow is NOT IMPLEMENTED YET."
        )
    else:
        raise ValueError(
            "WORKFLOW: the data type indicated in the sample sheet is neither DDA nor DIA."
        )
else:
    if path.exists(wf_path_in):
        local_path = wf_path_in
    else:
        raise FileNotFoundError(f"Supplied workfkow path '{wf_path_in}' is not valid.")


# import workflow and add path to database
with open(local_path, "r") as wf_file:
    wf = wf_file.read()
    wf = wf + f"\n{db_var}={path.abspath(db_path)}"

# export workflow
with open(wf_path_out, "w") as wf_out:
    wf_out.write(wf)

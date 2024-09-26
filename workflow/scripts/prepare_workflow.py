#!/usr/bin/python3

# FRAGPIPE WORKFLOWS based on sample type
# -----------------------------------------------------------------------------
#
# default workflow: LFQ-MBR (for DDA, data dependent acqisition)
#
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
default_wfs = snakemake.params["workflow_dir"]
db_path = snakemake.input["database"]
sample_sheet = snakemake.input["samplesheet"]
output_log = snakemake.log["path"]
db_var = "database.db-path"
local_path = ""
log = []
error = []

# determine if workflow file was supplied
if wf_path_in == "from_samplesheet":
    # if not, determine wf from sample type
    data_type = []
    with open(sample_sheet, "r") as sample_file:
        for sf in sample_file.readlines():
            data_type += [sf.split("\t")[3]]
    for dt in list(set(data_type)):
        if dt == "DDA":
            local_path = path.join(default_wfs, "LFQ-MBR.workflow")
            log += [f"Detected DDA samples, choosing default workflow: {local_path}"]
        elif dt == "DIA":
            error += ["Detected DIA samples, a default workflow is not implemented yet"]
        else:
            error += [
                "The data type indicated in the sample sheet is neither DDA nor DIA"
            ]
else:
    if path.exists(wf_path_in):
        local_path = wf_path_in
    else:
        error += [f"Supplied workfkow path '{wf_path_in}' is not a valid path"]

if path.exists(local_path):
    # import workflow and add path to database
    with open(local_path, "r") as wf_file:
        wf = wf_file.read()
        wf = wf + f"\n{db_var}={path.abspath(db_path)}"
    # export workflow
    with open(wf_path_out, "w") as wf_out:
        wf_out.write(wf)
    log += ["Added database entry to workflow"]

# print error/log messages
if error:
    print("\n".join(error))
    raise ValueError(
        "Location or format of the supplied workflow was not correct, quitting"
    )
else:
    log += [f"Module finished successfully"]
    log = ["WORKFLOW: " + i for i in log]
    with open(output_log, "w") as log_file:
        log_file.write("\n".join(log))

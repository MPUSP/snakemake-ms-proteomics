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

from pathlib import Path
import pandas as pd

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
    data_type = pd.read_csv(sample_sheet, header=None, sep="\t").iloc[:, 3].unique()
    for dt in list(set(data_type)):
        if dt == "DDA":
            local_path = Path(default_wfs) / "LFQ-MBR.workflow"
            log += [f"Detected DDA samples, choosing default workflow: {local_path}"]
        elif dt == "DIA":
            local_path = Path(default_wfs) / "DIA_SpecLib_Quant.workflow"
            log += [f"Detected DIA samples, choosing default workflow: {local_path}"]
        else:
            error += [
                "The data type indicated in the sample sheet is none of: 'DDA', 'DIA'"
            ]
else:
    if Path(wf_path_in).exists():
        local_path = wf_path_in
    else:
        error += [f"Supplied workflow path '{wf_path_in}' is not a valid path"]

if Path(local_path).exists():
    # import workflow and add path to database
    with open(local_path, "r") as wf_file:
        wf = wf_file.read()
        wf = wf + f"\n{db_var}={Path(db_path).resolve().as_posix()}\n"
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

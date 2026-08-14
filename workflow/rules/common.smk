# import modules
import pandas as pd
from snakemake.utils import validate
from pathlib import Path

# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)


# validate sample sheet and config file
validate(samples, schema="../../workflow/schemas/samples.schema.yml")
validate(config, schema="../../workflow/schemas/config.schema.yml")


def wfpath(file):
    wf = Path(workflow.basedir) / file
    return wf.resolve().as_posix()

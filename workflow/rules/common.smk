# import basic packages
from os import path
from os import listdir


# construct target paths
def out(file):
    outpath = path.join(config["output"]["path"], file)
    return path.abspath(outpath)


def wfpath(file):
    wf = path.join(workflow.basedir, file)
    return path.abspath(wf)


# print input parameters
def print_params(config):
    print("\n +++ WORKFLOW PARAMETERS +++ \n")
    for i in config.keys():
        print(f"  - {i}: {config.get(i)}")
    print("\n +++++++++++++++++++++++++++ \n")

# import basic packages
from os import path
from os import listdir

# construct target paths
def out(file):
    return path.join(config["output"]["path"], file)

def wfpath(file):
    return path.join(workflow.basedir, file)

# print input parameters
def print_params(config):
    print("\n +++ WORKFLOW PARAMETERS +++ \n")
    for i in config.keys():
        print(f"  - {i}: {config.get(i)}")
    print("\n +++++++++++++++++++++++++++ \n")

from pathlib import Path
import logging
OUTPUT_PATH = "ikarus/python/test/"

# @todo define a common place where python tests should output the output files to
def output_path():
    cwd = Path.cwd()

    if cwd.name == "test":
        return ""
    else:
        return OUTPUT_PATH


# SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from pathlib import Path

OUTPUT_PATH = "ikarus/python/test/"


# @todo define a common place where python tests should output the output files to
def output_path():
    cwd = Path.cwd()

    if cwd.name == "test":
        return ""
    else:
        return OUTPUT_PATH

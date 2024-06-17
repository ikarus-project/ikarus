# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
import os
os.environ["DUNE_LOG_LEVEL"] = "debug"
os.environ["DUNE_SAVE_BUILD"] = "terminal"
os.environ["DUNE_CMAKE_FLAGS"] = "-CMAKE_BUILD_TYPE=debug"
def setDebugFlags():
    from ikarus.generator import MySimpleGenerator
    MySimpleGenerator.setFlags("-g -Wfatal-errors",noChecks=True)

def unsetDebugFlags():
    from ikarus.generator import MySimpleGenerator
    MySimpleGenerator.reset()

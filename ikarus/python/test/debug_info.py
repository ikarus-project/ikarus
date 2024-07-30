# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
import os

os.environ["DUNE_LOG_LEVEL"] = "debug"
os.environ["DUNE_SAVE_BUILD"] = "terminal"


def setDebugFlags():
    import dune.generator as generator

    generator.setFlags("-g ", noChecks=False)


def unsetDebugFlags():
    import dune.generator as generator

    generator.reset()

# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
import os
import dune.generator as generator


os.environ["DUNE_LOG_LEVEL"] = "debug"
os.environ["DUNE_SAVE_BUILD"] = "terminal"


def setDebugFlags():
    def apply_debug_flags():
        generator.setFlags("-g ", noChecks=False)

    # Check if the environment variable is set
    build_type = os.environ.get("IKARUS_PYTHON_TEST_BUILD_TYPE")

    # Apply flags only if the build type is Debug, or if the variable is not set
    if build_type == "Debug" or build_type is None:
        apply_debug_flags()
        print("Python Generator: Debug Flag enabled")


def unsetDebugFlags():
    import dune.generator as generator

    generator.reset()

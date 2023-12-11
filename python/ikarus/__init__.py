# SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
try:
    from dune.packagemetadata import registerExternalModule
    import pathlib

    # register ikarus to be recognized by dune-py (code generation module)
    # as a module of the dune universe
    registerExternalModule(
        moduleName="ikarus",
        modulePath=str(pathlib.Path(__file__).parent.resolve()),
    )

except ImportError:
    pass

from ._ikarus import *
from .dirichlet_values import dirichletValues
from .basis import basis

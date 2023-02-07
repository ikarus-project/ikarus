# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
from ._toVoigt import *
__all__ = ["add","to_voigt"]




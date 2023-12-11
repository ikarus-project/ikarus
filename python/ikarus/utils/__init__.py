# SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator


def boundaryPatch(gridView, booleanVector):
    generator = MySimpleGenerator("BoundaryPatch", "Ikarus::Python")
    element_type = f"BoundaryPatch<{gridView.cppTypeName}>"

    includes = []
    includes += ["dune/fufem/boundarypatch.hh"]
    includes += ["ikarus/python/utils/boundarypatch.hh"]
    includes += gridView._includes
    moduleName = "boundaryPatch_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return module.BoundaryPatch(gridView, booleanVector)

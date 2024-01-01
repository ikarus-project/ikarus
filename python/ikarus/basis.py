# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from dune.functions.globalbasis import preBasisTypeName
from dune.functions import defaultGlobalBasis
from .generator import MySimpleGenerator


def basis(gv, tree):
    generator = MySimpleGenerator("Basis", "Ikarus::Python")

    pbfName = preBasisTypeName(tree, gv.cppTypeName)
    element_type = f"Ikarus::Basis<{pbfName}>"

    includes = []
    includes += ["ikarus/python/basis/basis.hh"]
    moduleName = "Basis_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    basis = defaultGlobalBasis(gv, tree)
    return module.Basis(basis)

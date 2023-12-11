# SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


from dune.common.hashit import hashIt
from .generator import MySimpleGenerator


def dirichletValues(basis):
    generator = MySimpleGenerator("DirichletValues", "Ikarus::Python")
    element_type = f"Ikarus::DirichletValues<{basis.cppTypeName},Eigen::VectorX<bool>>"

    includes = []
    includes += basis._includes
    includes += ["ikarus/assembler/simpleassemblers.hh"]
    includes += ["ikarus/python/dirichletvalues/dirichletvalues.hh"]
    moduleName = "SparseFlatAssembler_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return module.DirichletValues(basis)

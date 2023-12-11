# SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator


def sparseFlatAssembler(fes, dirichletValues):
    element_type = f"Ikarus::SparseFlatAssembler<std::vector<{fes[0].cppTypeName}>,{dirichletValues.cppTypeName}>"
    generator = MySimpleGenerator("SparseFlatAssembler", "Ikarus::Python")

    includes = []
    includes += ["ikarus/assembler/simpleassemblers.hh"]
    includes += fes[0]._includes  # include header of finite element
    includes += ["ikarus/python/assembler/flatassembler.hh"]
    moduleName = "SparseFlatAssembler_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return module.SparseFlatAssembler(fes, dirichletValues)


def denseFlatAssembler(fes, dirichletValues):
    element_type = f"Ikarus::DenseFlatAssembler<std::vector<{fes[0].cppTypeName}>,{dirichletValues.cppTypeName}>"
    generator = MySimpleGenerator("DenseFlatAssembler", "Ikarus::Python")

    includes = []
    includes += ["ikarus/assembler/simpleassemblers.hh"]
    includes += ["ikarus/finiteelements/mechanics/linearelastic.hh"]
    includes += ["ikarus/python/assembler/flatassembler.hh"]
    includes += fes[0]._includes  # include header of finite element
    moduleName = "SparseFlatAssembler_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return module.DenseFlatAssembler(fes, dirichletValues)

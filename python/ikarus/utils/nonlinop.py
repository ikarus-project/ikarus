# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator
from dune.generator.algorithm import cppType


def NonLinearOperator(functions, parameters):

    generator = MySimpleGenerator("NonLinearOperator", "Ikarus::Python")
    includes = []
    element_type = "Ikarus::NonLinearOperator<Ikarus::Impl::Functions<"
    # to create the correct function we have to deduce the return type of the given finctions
    from typing import get_type_hints

    for function in functions:
        value = function(*parameters)
        returnValueType, i = cppType(value)
        element_type += "const std::function<" + returnValueType + "("
        for par in parameters:
            t, i = cppType(par)
            element_type += "const " + t + ","
            includes += i
        element_type = element_type.strip(",") + ")>&,"

    element_type = element_type.strip(",") + ">,Ikarus::Impl::Parameter<"

    for par in parameters:
        t, i = cppType(par)
        element_type += t + ","
        includes += i
    element_type = element_type.strip(",") + ">>"

    includes += ["ikarus/python/utils/nonlinearoperator.hh"]
    moduleName = "NonLinearOperator_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return module.NonLinearOperator(functions, parameters)


def makeNonLinearOperator(assembler, requirement=None,affordances=None, enforcingBCOption=None):
    generator = MySimpleGenerator("NonLinearOperatorFactory", "Ikarus::Python")
    includes = []
    element_type = f"Ikarus::Python::NonLinearOperatorFactoryWrapper<std::shared_ptr<{assembler.cppTypeName}>>"

    includes += ["ikarus/python/utils/nonlinearoperator.hh"]
    includes+= assembler.cppIncludes
    moduleName = "NonLinearOperatorFactory_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    factory= module.NonLinearOperatorFactory(assembler)
    factory
    element_type = factory.nonLinearOperatorType
    moduleNameNonLinOp = "NonLinearOperator" + hashIt(element_type)
    moduleForNonLinOp = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return factory.op(requirement,affordances,enforcingBCOption)

# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


from io import StringIO
import types
from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator
from dune.generator.algorithm import cppType, run

def __SubOperator(nonLinOp):
    # def __subOperatorFunc(nonLinOp,*args: int):
    #     runCode="""
    #     template <class NLO>
    #     void subOperator(const NLO& nlo)
    #     {{
    #     return nlo.template subOperator<{indices}>();
    #     }}
    #     """.format(indices=",".join([str(i) for i in args]))
    #     return run("subOperator",StringIO(runCode),nonLinOp)

    def __subOperatorFunc(nonLinOp, *args):
        sorted = all(args[i] <= args[i+1] for i in range(len(args) - 1))
        if(not sorted):
            raise ValueError("The indices passed to the subOperator function are not sorted")
        if(len(args)>=nonLinOp.numberOfFunctions):
            raise ValueError("You passed too many arguments to the subOperator function")
        elif(len(args)==0):
            raise ValueError("At least one argument is needed")
        elif(any([arg>=nonLinOp.numberOfFunctions-1 for arg in args])):
            raise ValueError("You passed an index that is too high")
        elif(len(args)==1):
            if(args[0] == 0):
                return nonLinOp.__subOperator0()
            elif(args[0] == 1):
                return nonLinOp.__subOperator1()
            elif(args[0] == 2):
                return nonLinOp.__subOperator2()
        elif(len(args)==2):
            if(args[0] == 0 and args[1] == 1):
                return nonLinOp.__subOperator01()
            elif(args[0] == 1 and args[1] == 2):
                return nonLinOp.__subOperator12()
        else:
            raise ValueError("The subOperator function does not know how to handle the given indices")
    return __subOperatorFunc


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
        includes=includes, typeName=element_type, moduleName=moduleName,dynamicAttr=True # dynamicAttr is needed to allow the addition of the subOperator function
    )
    nonLinOp = module.NonLinearOperator(functions, parameters)

    nonLinOp.subOperator = types.MethodType(__SubOperator(nonLinOp),nonLinOp)
    return nonLinOp


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

    element_type = factory.nonLinearOperatorType
    # moduleNameNonLinOp = "NonLinearOperator" + hashIt(element_type)
    # moduleForNonLinOp = generator.load(
    #     includes=includes, typeName=element_type, moduleName=moduleName
    # )

    nonLinOp =  factory.op(requirement,affordances,enforcingBCOption)
    nonLinOp.subOperator = types.MethodType(__SubOperator(nonLinOp),nonLinOp)
    return nonLinOp
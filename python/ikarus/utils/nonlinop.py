# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


from io import StringIO
import types
from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator
from dune.generator.algorithm import cppType, run

def __SubOperator(nonLinOp):
    # This would be nice to have, but it is not possible to generate the code for the subOperator function
    def __subOperatorFunc(nonLinOp,*args: int):
        runCode="""
        template <class NLO>
        void subOperator( NLO& nlo)
        {{
        return nlo.template subOperator<{indices}>();
        }}
        """.format(indices=",".join([str(i) for i in args]))
        return run("subOperator",StringIO(runCode),nonLinOp)

    def __subOperatorFunc(nonLinOp, *args):
        from itertools import accumulate
        sorted = all(args[i] <= args[i+1] for i in range(len(args) - 1))
        if(not sorted):
            raise ValueError("The indices passed to the subOperator function are not sorted")
        if(len(args)>=nonLinOp.numberOfFunctions):
            raise ValueError("You passed too many arguments to the subOperator function")
        elif(len(args)==0):
            raise ValueError("At least one argument is needed")
        elif(any([arg>nonLinOp.numberOfFunctions-1 for arg in args])):
            raise ValueError("You passed an index that is too high")

        func = getattr(nonLinOp, '__subOperator'.join([str(x) for x in args]))
        subOp = func()
        subOp.subOperator = types.MethodType(__SubOperator(subOp),subOp)
        return subOp
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


def __op(nonLinOpFactory):
    # This would be nice to have, but it is not possible to generate the code for the subOperator function
    def __opFunc(nonLinOpFactory,requirement,affordances,enforcingBCOption,*args: int):
        runCode="""
        template <class NLO>
        void subOperator( NLO& nlo)
        {{
        return nlo.template subOperator<{indices}>();
        }}
        """.format(indices=",".join([str(i) for i in args]))
        op = run("subOperator",StringIO(runCode),nonLinOp)
        op.subOperator = types.MethodType(__opFunc(op),op)
        return op

    # def __opFunc(nonLinOpFactory, *args):
    #     from itertools import accumulate
    #     sorted = all(args[i] <= args[i+1] for i in range(len(args) - 1))
    #     if(not sorted):
    #         raise ValueError("The indices passed to the subOperator function are not sorted")
    #     if(len(args)>=nonLinOp.numberOfFunctions):
    #         raise ValueError("You passed too many arguments to the subOperator function")
    #     elif(len(args)==0):
    #         raise ValueError("At least one argument is needed")
    #     elif(any([arg>nonLinOp.numberOfFunctions-1 for arg in args])):
    #         raise ValueError("You passed an index that is too high")

    return __opFunc

def makeNonLinearOperator(assembler,derivativeIndices=None, requirement=None,affordances=None, enforcingBCOption=None,*args: int):
    generator = MySimpleGenerator("NonLinearOperatorFactory", "Ikarus::Python")
    includes = []
    element_type = f"Ikarus::Python::NonLinearOperatorFactoryWrapper<std::shared_ptr<{assembler.cppTypeName}>>"

    includes += ["ikarus/python/utils/nonlinearoperator.hh"]
    includes+= assembler.cppIncludes
    moduleName = "NonLinearOperatorFactory_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName, dynamicAttr=True
    )
    factory= module.NonLinearOperatorFactory(assembler)
    factory.op = types.MethodType(__op(factory),factory)
    nonLinOp =  factory.op(requirement,affordances,enforcingBCOption,*args)
    nonLinOp.subOperator = types.MethodType(__SubOperator(nonLinOp),nonLinOp)
    return nonLinOp
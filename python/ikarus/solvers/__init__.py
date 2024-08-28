# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator
from typing import Callable


def TrustRegion(nonlinearOperator,PreConditioner="IncompleteCholesky", updateFunction=None):
    """
    @brief 

    @param 

    @return: 
    """
    generator = MySimpleGenerator("TrustRegion", "Ikarus::Python")
    #includes = ["ikarus/solver/nonlinearsolver/trustregion.hh"]
    includes = ["ikarus/python/solvers/registertrustregion.hh"]
    includes += nonlinearOperator.cppIncludes
    updateFunctionType = f"std::function<void({nonlinearOperator.firstParameterCppTypeName}&,const {nonlinearOperator.derivativeCppTypeName}&)>"
    element_type = (
        f"Ikarus::TrustRegion<{nonlinearOperator.cppTypeName},Ikarus::PreConditioner::{PreConditioner},{updateFunctionType}>"
    )
    print("element_type: ", element_type)

    moduleName = "TrustRegion" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    if updateFunction is None:
        return module.TrustRegion(nonlinearOperator)
    else:
        return module.TrustRegion(nonlinearOperator, updateFunction)
    

def NewtonRaphson(nonlinearOperator,linearSolver="UmfPack", updateFunction=None):
    """
    @brief 

    @param 

    @return: 
    """
    generator = MySimpleGenerator("NewtonRaphson", "Ikarus::Python")
    #includes = ["ikarus/solver/nonlinearsolver/trustregion.hh"]
    includes = ["ikarus/python/solvers/registernewton.hh"]
    includes += nonlinearOperator.cppIncludes
    if nonlinearOperator.numberOfFunctions == 2:
        updateFunctionType = f"std::function<void({nonlinearOperator.firstParameterCppTypeName}&,const {nonlinearOperator.valueCppTypeName}&)>"
    elif nonlinearOperator.numberOfFunctions == 3:
        updateFunctionType = f"std::function<void({nonlinearOperator.secondParameterCppTypeName}&,const {nonlinearOperator.derivativeCppTypeName}&)>"

    element_type = (
        f"Ikarus::NewtonRaphson<{nonlinearOperator.cppTypeName},Ikarus::LinearSolver,{updateFunctionType}>"
    )
    print("element_type: ", element_type)

    moduleName = "NewtonRaphson" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    if updateFunction is None:
        return module.NewtonRaphson(nonlinearOperator,linearSolver)
    else:
        return module.NewtonRaphson(nonlinearOperator,linearSolver, updateFunction)

# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator
from typing import Callable
#from . import PreConditioner

def __updateFunctionType(nonlinearOperator):
    assert "Eigen::Ref" in nonlinearOperator.firstParameterCppTypeName
    print("nonlinearOperator.numberOfFunctions: ", nonlinearOperator.numberOfFunctions)
    if nonlinearOperator.numberOfFunctions == 2:
        updateFunctionType = f"std::function<void({nonlinearOperator.firstParameterCppTypeName},const Eigen::Ref<const {nonlinearOperator.valueCppTypeName}>&)>"
    elif nonlinearOperator.numberOfFunctions == 3:
        updateFunctionType = f"std::function<void({nonlinearOperator.secondParameterCppTypeName},const Eigen::Ref<const {nonlinearOperator.derivativeCppTypeName}>&)>"
    return updateFunctionType


def TrustRegion(
    nonlinearOperator,preConditioner=None, updateFunction=None
):
    """
      @brief Constructs and returns a TrustRegion solver objects, providing an interface to the C++ `TrustRegion` class
             using Pybind11.

      This solver applies the Trust Region algorithm, typically used for solving nonlinear optimization problems,
      with options for customization through preconditioning and update functions.

    @param nonlinearOperator A nonlinear operator instance that defines the problem to be solved.
    This object should have a method to evaluate the cost function and its first and second derivatives.

      @param PreConditioner Specifies the preconditioner type to be used in the solver. Default is `solvers.PreConditioner.IncompleteCholesky`.
             This string should match a valid preconditioner in the `Ikarus::PreConditioner` enum. Optional.

      @param updateFunction A custom update function for the TrustRegion solver. It should accept two parameters:
             - The first parameter with the same type as `nonlinearOperator.firstParameterCppTypeName`
             - The second parameter of type `Eigen::Ref<const {nonlinearOperator.derivativeCppTypeName}>`
             Optional.

      @return A Pybind11-wrapped TrustRegion solver instance. If `updateFunction` is provided, the solver will
              incorporate the specified update function. If the `updateFunction` is `None`, the solver will use a default configuration with an update function. See on the c++ side for details

      @exception AssertionError If `nonlinearOperator.firstParameterCppTypeName` does not contain `"Eigen::Ref"`,
                 indicating an incompatibility in the expected parameter type.
    """
    generator = MySimpleGenerator("TrustRegion", "Ikarus::Python")
    # includes = ["ikarus/solver/nonlinearsolver/trustregion.hh"]
    includes = ["ikarus/python/solvers/registertrustregion.hh"]
    includes += nonlinearOperator.cppIncludes
    assert "Eigen::Ref" in nonlinearOperator.firstParameterCppTypeName
    assert nonlinearOperator.numberOfFunctions==3

    updateFunctionType = f"std::function<void( {nonlinearOperator.firstParameterCppTypeName},const Eigen::Ref<const {nonlinearOperator.derivativeCppTypeName}>&)>"

    if preConditioner is None:
        preConditioner = PreConditioner.IncompleteCholesky
    element_type = f"Ikarus::TrustRegion<{nonlinearOperator.cppTypeName},Ikarus::PreConditioner::{str(preConditioner).replace('PreConditioner.', '')},{updateFunctionType}>"

    moduleName = "TrustRegion" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    if updateFunction is None:
        return module.TrustRegion(nonlinearOperator)
    else:
        return module.TrustRegion(nonlinearOperator, updateFunction)





def NewtonRaphson(nonlinearOperator, linearSolver=None, updateFunction=None):
    """
     @brief Constructs and returns a NewtonRaphson solver objects, providing an interface to the C++ `NewtonRaphson` class
             using Pybind11.

      This solver applies the Trust Region algorithm, typically used for solving nonlinear optimization problems,
      with options for customization through preconditioning and update functions.

    @param nonlinearOperator A nonlinear operator instance that defines the problem to be solved.
    This object should have a method to evaluate the residial function and its jacobian.

      @param linearSolver Specifies the linear solver used.

      @param updateFunction A custom update function for the TrustRegion solver. It should accept two parameters:
             - The first parameter with the same type as `nonlinearOperator.firstParameterCppTypeName`
             - The second parameter of type `Eigen::Ref<const {nonlinearOperator.derivativeCppTypeName}>`
             Optional.

      @return A Pybind11-wrapped NewtonRaphson solver instance. If `updateFunction` is provided, the solver will
              incorporate the specified update function. If the `updateFunction` is `None`, the solver will use a default configuration with an update function. See on the c++ side for details

    """
    generator = MySimpleGenerator("NewtonRaphson", "Ikarus::Python")
    includes = ["ikarus/python/solvers/registernewton.hh"]
    includes += nonlinearOperator.cppIncludes
    updateFunctionType = __updateFunctionType(nonlinearOperator)

    element_type = f"Ikarus::NewtonRaphson<{nonlinearOperator.cppTypeName},Ikarus::LinearSolver,{updateFunctionType}>"
    print("element_type: ", element_type)

    moduleName = "NewtonRaphson" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    if updateFunction is None:
        return module.NewtonRaphson(nonlinearOperator, linearSolver)
    else:
        return module.NewtonRaphson(nonlinearOperator, linearSolver, updateFunction)


def NewtonRaphsonWithSubsidiaryFunction(
    nonlinearOperator, linearSolver=None, updateFunction=None
):
    """
     @brief Constructs and returns a NewtonRaphsonWithSubsidiaryFunction solver objects, providing an interface to the C++ `NewtonRaphsonWithSubsidiaryFunction` class
             using Pybind11.

      This solver applies the Trust Region algorithm, typically used for solving nonlinear optimization problems,
      with options for customization through preconditioning and update functions.

    @param nonlinearOperator A nonlinear operator instance that defines the problem to be solved.
    This object should have a method to evaluate the residial function and its jacobian.

      @param linearSolver Specifies the linear solver used.

      @param updateFunction A custom update function for the TrustRegion solver. It should accept two parameters:
             - The first parameter with the same type as `nonlinearOperator.firstParameterCppTypeName`
             - The second parameter of type `Eigen::Ref<const {nonlinearOperator.derivativeCppTypeName}>`
             Optional.

      @return A Pybind11-wrapped NewtonRaphsonWithSubsidiaryFunction solver instance. If `updateFunction` is provided, the solver will
              incorporate the specified update function. If the `updateFunction` is `None`, the solver will use a default configuration with an update function. See on the c++ side for details

    """
    generator = MySimpleGenerator(
        "NewtonRaphsonWithSubsidiaryFunction", "Ikarus::Python"
    )
    includes = ["ikarus/python/solvers/registernewtonwithsubsidaryfunction.hh"]
    includes += nonlinearOperator.cppIncludes
    updateFunctionType = __updateFunctionType(nonlinearOperator)

    element_type = f"Ikarus::NewtonRaphsonWithSubsidiaryFunction<{nonlinearOperator.cppTypeName},Ikarus::LinearSolver,{updateFunctionType}>"
    print("element_type: ", element_type)

    moduleName = "NewtonRaphsonWithSubsidiaryFunction" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    if updateFunction is None:
        return module.NewtonRaphsonWithSubsidiaryFunction(
            nonlinearOperator, linearSolver
        )
    else:
        return module.NewtonRaphsonWithSubsidiaryFunction(
            nonlinearOperator, linearSolver, updateFunction
        )

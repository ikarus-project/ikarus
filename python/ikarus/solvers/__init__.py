# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator
from typing import Callable


def TrustRegion(nonlinearOperator,PreConditioner="IncompleteCholesky", updateFunction=None):
    """
    @brief Registers a pre-element with the specified parameters.

    @param name: The name of the pre-element.
    @param includes: List of additional include files.
    @param element_type: The type of the pre-element.
    @param args: Additional arguments required for the pre-element registration.

    @return: The registered pre-element function.
    """
    generator = MySimpleGenerator("TrustRegion", "Ikarus::Python")
    #includes = ["ikarus/solver/nonlinearsolver/trustregion.hh"]
    includes = ["ikarus/python/solvers/registersolver.hh"]
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

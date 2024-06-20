# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from .generator import MySimpleGenerator

def dirichletValues(basis):
    """
    @brief Creates a Dirichlet values handler for the given basis.

    @param basis: The basis.

    @return: The created Dirichlet values handler.
    """
    generator = MySimpleGenerator("DirichletValues", "Ikarus::Python")
    element_type = f"Ikarus::DirichletValues<{basis.cppTypeName},Eigen::VectorX<bool>>"

    includes = []
    includes += basis._includes
    includes += ["ikarus/utils/dirichletvalues.hh"]
    includes += ["ikarus/python/dirichletvalues/dirichletvalues.hh"]
    moduleName = "dirichletValues_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return module.DirichletValues(basis)

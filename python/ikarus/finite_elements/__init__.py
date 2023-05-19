# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator


def linearElasticElement(
    basis, element, youngsMod, nu, volumeLoad=None, bp=None, neumannBoundaryLoad=None
):
    if not (
        (bp is None and neumannBoundaryLoad is None)
        or (bp is not None and neumannBoundaryLoad is not None)
    ):
        raise TypeError(
            "If you provide a boundary patch you should also provide a boundary load!"
        )

    generator = MySimpleGenerator("LinearElastic", "Ikarus::Python")
    element_type = f"Ikarus::LinearElastic<{basis.cppTypeName},Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>"

    includes = []
    includes += ["ikarus/finiteElements/mechanics/linearElastic.hh"]
    includes += ["ikarus/python/finiteElements/linearElastic.hh"]
    includes += basis._includes
    includes += element._includes
    moduleName = "linearElastic_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    # https://pybind11.readthedocs.io/en/stable/advanced/functions.html#allow-prohibiting-none-arguments
    if volumeLoad is None:
        return module.LinearElastic(basis, element, youngsMod, nu)
    elif bp is None and neumannBoundaryLoad is None:
        return module.LinearElastic(basis, element, youngsMod, nu, volumeLoad)
    else:
        return module.LinearElastic(
            basis, element, youngsMod, nu, volumeLoad, bp, neumannBoundaryLoad
        )


def nonLinearElasticElement(
    basis, element, material, volumeLoad=None, bp=None, neumannBoundaryLoad=None
):
    if not (
        (bp is None and neumannBoundaryLoad is None)
        or (bp is not None and neumannBoundaryLoad is not None)
    ):
        raise TypeError(
            "If you provide a boundary patch you should also provide a boundary load!"
        )

    generator = MySimpleGenerator("NonLinearElastic", "Ikarus::Python")
    element_type = f"Ikarus::NonLinearElastic<{basis.cppTypeName},  {material.cppTypeName} ,Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>"

    includes = []
    includes += ["ikarus/finiteElements/mechanics/nonLinearElastic.hh"]
    includes += ["ikarus/python/finiteElements/nonLinearElastic.hh"]
    includes += basis._includes
    includes += material._includes
    includes += element._includes
    moduleName = "nonLinearElastic_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    # https://pybind11.readthedocs.io/en/stable/advanced/functions.html#allow-prohibiting-none-arguments
    if volumeLoad is None:
        return module.NonLinearElastic(basis, element, material)
    elif bp is None and neumannBoundaryLoad is None:
        return module.NonLinearElastic(basis, element, material, volumeLoad)
    else:
        return module.NonLinearElastic(
            basis, element, material, volumeLoad, bp, neumannBoundaryLoad
        )

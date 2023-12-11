# SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator


def elementConstructorDecoratorFactory(
    allowsMaterial, elementHeader, autodiffWrapper=None
):
    def decorator(func):
        def wrapper(
            basis,
            element,
            youngsMod,
            nu,
            volumeLoad=None,
            bp=None,
            neumannBoundaryLoad=None,
        ):
            generator = MySimpleGenerator(func.__name__, "Ikarus::Python")
            # if allowsMaterial:
            includes = [elementHeader]
            if autodiffWrapper is None:
                element_type = (
                    "Ikarus::"
                    + func.__name__
                    + f"<{basis.cppTypeName},Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>"
                )
            else:
                includes += ["ikarus/finiteelements/febases/autodifffe.hh"]
                element_type = (
                    "Ikarus::"
                    + autodiffWrapper
                    + "<Ikarus::"
                    + func.__name__
                    + f"<{basis.cppTypeName},Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>,"
                    f"Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>"
                )

            # else:
            #     element_type = "Ikarus::" + func.__name__ + f"<{basis.cppTypeName},  {material.cppTypeName} ,Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>"

            if not (
                (bp is None and neumannBoundaryLoad is None)
                or (bp is not None and neumannBoundaryLoad is not None)
            ):
                raise TypeError(
                    "If you provide a boundary patch you should also provide a boundary load!"
                )

            includes += basis._includes
            # includes += material._includes
            includes += element._includes
            moduleName = func.__name__ + "_" + hashIt(element_type)
            module = generator.load(
                includes=includes, typeName=element_type, moduleName=moduleName
            )
            # https://pybind11.readthedocs.io/en/stable/advanced/functions.html#allow-prohibiting-none-arguments
            if volumeLoad is None:
                return eval(
                    "module." + func.__name__ + "(basis, element, youngsMod, nu)"
                )
            elif bp is None and neumannBoundaryLoad is None:
                return eval(
                    "module."
                    + func.__name__
                    + "(basis, element, youngsMod, nu, volumeLoad)"
                )
            else:
                return eval(
                    "module."
                    + func.__name__
                    + "(basis, element, youngsMod, nu, volumeLoad, bp, neumannBoundaryLoad)"
                )

        return wrapper

    return decorator


@elementConstructorDecoratorFactory(
    False, "ikarus/python/finiteelements/linearelastic.hh"
)
def LinearElastic(
    basis, element, youngsMod, nu, volumeLoad=None, bp=None, neumannBoundaryLoad=None
):
    return elementConstructorDecoratorFactory(
        False, "ikarus/python/finiteelements/linearelastic.hh"
    )


def KirchhoffLoveShell(
    basis,
    element,
    youngsMod,
    nu,
    thickness,
    volumeLoad=None,
    bp=None,
    neumannBoundaryLoad=None,
):
    generator = MySimpleGenerator("KirchhoffLoveShell", "Ikarus::Python")
    includes = ["ikarus/python/finiteelements/linearelastic.hh"]

    includes += ["ikarus/finiteelements/febases/autodifffe.hh"]
    autodiffWrapper = "AutoDiffFE"
    element_type = (
        "Ikarus::"
        + autodiffWrapper
        + f"<Ikarus::KirchhoffLoveShell<{basis.cppTypeName},Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>,"
        f"Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>"
    )

    # else:
    #     element_type = "Ikarus::" + func.__name__ + f"<{basis.cppTypeName},  {material.cppTypeName} ,Ikarus::FErequirements<Eigen::Ref<Eigen::VectorXd>>,true>"

    if not (
        (bp is None and neumannBoundaryLoad is None)
        or (bp is not None and neumannBoundaryLoad is not None)
    ):
        raise TypeError(
            "If you provide a boundary patch you should also provide a boundary load!"
        )

    includes += basis._includes
    # includes += material._includes
    includes += element._includes
    moduleName = "KirchhoffLoveShell" + "_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    # https://pybind11.readthedocs.io/en/stable/advanced/functions.html#allow-prohibiting-none-arguments
    if volumeLoad is None:
        return module.KirchhoffLoveShell(basis, element, youngsMod, nu, thickness)
    elif bp is None and neumannBoundaryLoad is None:
        return module.KirchhoffLoveShell(
            basis, element, youngsMod, nu, thickness, volumeLoad
        )
    else:
        return module.KirchhoffLoveShell(
            basis,
            element,
            youngsMod,
            nu,
            thickness,
            volumeLoad,
            bp,
            neumannBoundaryLoad,
        )


def NonLinearElastic(
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
    includes += ["ikarus/python/finiteelements/nonlinearelastic.hh"]
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

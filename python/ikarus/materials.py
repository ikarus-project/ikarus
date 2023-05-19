# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


from dune.common.hashit import hashIt
from .generator import MySimpleGenerator


def materialConstructorDecorator(func):
    def wrapper(emodul, nu):
        generator = MySimpleGenerator(func.__name__, "Ikarus::Python")
        element_type = "Ikarus::" + func.__name__ + "<double>"

        includes = []
        includes += ["ikarus/finiteElements/mechanics/materials.hh"]
        includes += ["ikarus/python/finiteElements/materials/material.hh"]
        moduleName = func.__name__ + "_" + hashIt(element_type)
        module = generator.load(
            includes=includes, typeName=element_type, moduleName=moduleName
        )

        return eval("module." + func.__name__ + "(emodul,nu)")

    return wrapper


@materialConstructorDecorator
def StVenantKirchhoff(emodul, nu):
    return materialConstructorDecorator


@materialConstructorDecorator
def NeoHooke(emodul, nu):
    return materialConstructorDecorator

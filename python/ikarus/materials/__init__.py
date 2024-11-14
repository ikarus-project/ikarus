# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt

from ikarus.generator import MySimpleGenerator
from dune.vtk import FormatTypes, DataTypes
import dune.vtk

import types
from enum import Enum
from typing import Union

# The list of supported dataCollectors
MuesliElastic = Enum("MuesliElastic", ["LinearElasticity"])
MuesliFinite = Enum(
    "MuesliFinite",
    [
        "StVenantKirchhoff",
        "NeoHooke",
        "MooneyRivlin",
        "Yeoh",
        "ArrudaBoyce",
    ],
)


def muesliMaterial(materialTag: Union[MuesliElastic, MuesliFinite], **kwargs):
    assert  isinstance(materialTag, MuesliElastic) or isinstance(materialTag, MuesliFinite)

    materialTagSplitted = str(materialTag).split(".")
    wrappertype = materialTagSplitted[0]
    materialtype = f"Ikarus::Materials::Muesli::{materialTagSplitted[1]}"

    includes = []
    includes += ["ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh"]
    includes += ["ikarus/python/finiteelements/material.hh"]

    generator = MySimpleGenerator("MuesliMaterial", "Ikarus::Python")
    element_type = f"Ikarus::Materials::{wrappertype}<{materialtype}>"

    moduleName = "MuesliMaterial_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName,
    )

    return module.MuesliMaterial(**kwargs)

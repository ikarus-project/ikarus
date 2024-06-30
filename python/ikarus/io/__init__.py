# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator
from dune.vtk import FormatTypes, DataTypes

from logging import warn

dataCollectors = ["lagrange", "discontinous", "iga"]


def vtkWriter(
    assembler,
    dataCollector: str = None,
    order: int = 1,
    structured=False,
    format=FormatTypes.binary,
    datatype=DataTypes.Float32,
    headertype=DataTypes.UInt32,
):
    """ """

    includes = []

    # If we have a structured gridwriter, we default to the yaspdatacollector on the c++ side
    if dataCollector is not None and structured:
        warn("DataCollector Argument ignored, as Writer is set to structured")
        dataCollector = None

    structuredStr = "true" if structured else "false"
    gridViewName = assembler.grid().cppTypeName
    dataCollectorName: str = ""
    
    if dataCollector is not None:
        if dataCollector == dataCollectors[0]:
            dataCollectorName = (
                f"Dune::Vtk::LagrangeDataCollector<{gridViewName}, {order}>"
            )
            includes += ["dune/vtk/datacollectors/lagrangedatacollector.hh"]
        elif dataCollector == dataCollectors[1]:
            dataCollectorName = f"Dune::Vtk::DiscontinuousDataCollector<{gridViewName}>"
            includes += ["dune/vtk/datacollectors/discontinuousdatacollector.hh"]
        elif dataCollector == dataCollectors[2]:
            dataCollectorName = f"Dune::Vtk::DiscontinuousIgaDataCollector<{gridViewName}>"
            includes += ["dune/iga/io/igadatacollector.hh"]
        else:
            warn(
                f"DataCollector Argument ignored, as DataCollector {dataCollector} is unkwown"
            )
            dataCollector = None

    generator = MySimpleGenerator("VtkWriter", "Ikarus::Python")
    if dataCollector is None:
        element_type = f"Ikarus::Vtk::Writer<{assembler.cppTypeName}, {structuredStr}>"
    else:
        element_type = f"Ikarus::Vtk::Writer<{assembler.cppTypeName}, {structuredStr}, {dataCollectorName}>"

    includes += ["ikarus/assembler/simpleassemblers.hh"]
    includes += ["ikarus/io/vtkwriter.hh"]
    includes += ["ikarus/python/io/vtkwriter.hh"]
    includes += assembler._includes
    moduleName = "vtkWriter_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    return module.VtkWriter(assembler, format, datatype, headertype)

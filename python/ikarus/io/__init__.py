# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator
from dune.vtk import FormatTypes, DataTypes, vtkWriter
import dune.vtk

from warnings import warn
from io import StringIO
from dune.generator.algorithm import run
import types

# The list of supported dataCollectors
dataCollectors = ["lagrange", "discontinuous", "iga"]


def __addInterpolation(writer, flag: int):
    typeStr = "Ikarus::Vtk::asPointData" if flag == 1 else "Ikarus::Vtk::asCellData"

    def __addInterpolationFunc(writer, vals_, basis, name: str, size: int):
        containerStr = f"Dune::FieldVector<double, {size}>" if size > 1 else "double"

        runCode = """
            #define EIGEN_DEFAULT_TO_ROW_MAJOR 1
            #include <ikarus/python/io/vtkwriter.hh>
            #include <dune/python/pybind11/eigen.h>
            template <typename Writer, typename Basis>
            void addInterpolation(Writer& writer, const auto& vals_, const Basis& basis, std::string name) {{
                auto vals = vals_.template cast<Eigen::VectorX<double>>();
                writer.template addInterpolation<Basis, {containerStr}>(std::move(vals), basis, name, {typeStr});
            }}  
        """.format(
            containerStr=containerStr, typeStr=typeStr
        )

        return run("addInterpolation", StringIO(runCode), writer, vals_, basis, name)

    return __addInterpolationFunc


def vtkWriter(
    assembler,
    dataCollector: str = None,
    order: int = 1,
    format=FormatTypes.binary,
    datatype=DataTypes.Float32,
    headertype=DataTypes.UInt32,
):
    includes = []
    includes += ["dune/python/vtk/writer.hh"]
    includes += ["dune/vtk/writers/unstructuredgridwriter.hh"]
    includes += [
        "dune/vtk/writers/unstructuredgridwriter.hh",
        "dune/vtk/datacollectors/lagrangedatacollector.hh",
    ]
    includes += ["ikarus/io/vtkwriter.hh"]
    includes += ["ikarus/assembler/simpleassemblers.hh"]
    includes += assembler._includes

    gridViewName = assembler.grid.cppTypeName
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
            dataCollectorName = (
                f"Dune::Vtk::DiscontinuousIgaDataCollector<{gridViewName}>"
            )
            includes += ["dune/iga/io/igadatacollector.hh"]
        else:
            warn(
                f"DataCollector Argument ignored, as DataCollector {dataCollector} is unknown"
            )
            dataCollector = None

    generator = MySimpleGenerator("VtkWriter", "Ikarus::Python")

    defaultManager = f"Ikarus::Vtk::DefaultVTKWriterManager<{gridViewName}>"
    vtkWriterName = ""
    element_type = ""

    if dataCollector is None:
        vtkWriterName = f"typename {defaultManager}::template DefaultVTKWriter<>"
        dataCollectorName = f"typename {defaultManager}::DefaultDataCollector"
        element_type = f"Ikarus::Vtk::Writer<{assembler.cppTypeName}, {dataCollectorName}, {vtkWriterName}>"

    else:
        vtkWriterName = (
            f"typename {defaultManager}::template DefaultVTKWriter<{dataCollectorName}>"
        )
        element_type = f"Ikarus::Vtk::Writer<{assembler.cppTypeName}, {dataCollectorName}, {vtkWriterName}>"

    # Register VTKWriter
    dune.vtk.load(includes, vtkWriterName)

    includes += ["ikarus/python/io/vtkwriter.hh"]
    moduleName = "vtkWriter_" + hashIt(element_type)
    module = generator.load(
        includes=includes,
        typeName=element_type,
        moduleName=moduleName,
        baseClasses=[vtkWriterName],
        dynamicAttr=True,
    )
    writerModule = module.VtkWriter(assembler, format, datatype, headertype)
    writerModule.addInterpolationAsPointData = types.MethodType(
        __addInterpolation(writerModule, 1), writerModule
    )
    writerModule.addInterpolationAsCellData = types.MethodType(
        __addInterpolation(writerModule, 2), writerModule
    )

    return writerModule

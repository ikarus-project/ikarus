# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from dune.generator.algorithm import run

from ikarus.generator import MySimpleGenerator
from dune.vtk import FormatTypes, DataTypes
import dune.vtk

from warnings import warn
from io import StringIO
import types
from enum import Enum

from ._io import *

# The list of supported dataCollectors
DataCollector = Enum("DataCollector", ["default", "lagrange", "discontinuous", "iga"])


def __addInterpolation(writer):
    def __addInterpolationFunc(writer, vals_, basis, name: str, dataTag):
        runCode = """
            #include <ikarus/python/io/vtkwriter.hh>
            #include <dune/python/pybind11/eigen.h>
            template <typename Writer, typename Basis>
            void addInterpolation(Writer& writer, const auto& vals_, const Basis& basis, std::string name, int dataTag) {{
                auto vals = vals_.template cast<Eigen::VectorX<double>>();
                writer.template addInterpolation(std::move(vals), basis, name,  static_cast<Ikarus::Vtk::DataTag>(dataTag));
            }}  
        """
        return run(
            "addInterpolation",
            StringIO(runCode),
            writer,
            vals_,
            basis,
            name,
            dataTag.value,
        )

    return __addInterpolationFunc


def vtkWriter(
    assembler,
    dataCollector: DataCollector = DataCollector.default,
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

    gridViewName = assembler.gridView.cppTypeName
    dataCollectorName: str = ""

    if dataCollector is not DataCollector.default:
        if dataCollector == DataCollector.lagrange:
            dataCollectorName = (
                f"Dune::Vtk::LagrangeDataCollector<{gridViewName}, {order}>"
            )
        elif dataCollector == DataCollector.discontinuous:
            dataCollectorName = f"Dune::Vtk::DiscontinuousDataCollector<{gridViewName}>"
            includes += ["dune/vtk/datacollectors/discontinuousdatacollector.hh"]
        elif dataCollector == DataCollector.iga:
            dataCollectorName = (
                f"Dune::Vtk::DiscontinuousIgaDataCollector<{gridViewName}>"
            )
            includes += ["dune/iga/io/igadatacollector.hh"]

    generator = MySimpleGenerator("VtkWriter", "Ikarus::Python")

    defaultManager = f"Ikarus::Vtk::DefaultVTKWriterManager<{gridViewName}>"
    vtkWriterName = ""
    element_type = ""

    if dataCollector is DataCollector.default:
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
    writerModule.addInterpolation = types.MethodType(
        __addInterpolation(writerModule), writerModule
    )

    return writerModule

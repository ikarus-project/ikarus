// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vtkwriter.hh
 * \brief Python bindings for VtkWriter
 */

#pragma once

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/stl_bind.h>
#include <dune/vtk/vtkwriter.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/io/vtkdatatag.hh>
#include <ikarus/io/vtkwriter.hh>

namespace Ikarus::Python {

/**
 * \brief Register Python bindings for a VtkWriter class. \n
 *
 * The registered VtkWriter class provides functionalities for writing VTK files from assembled data.
 * This class supports adding result data as cell or point data and configuring VTK output formats.
 *
 * This function registers the following methods for the VtkWriter class:
 * - `setFormat(type: dune.vtk.FormatTypes)`
 * - `setDatatype(type: dune.vtk.DataTypes)`
 * - `setHeadertype(type: dune.vtk.DataTypes)`
 * - `addAllResultsAsCellData()`
 * - `addAllResultsAsPointData()`
 * - `addResultAsCellData(resType: str)`
 * - `addResultAsPointData(resType: str)`
 * - `write(fileName)`
 * - `addInterPolation(writer, vals_::np.array, basis, name: str, size: int)`
 * - `addPointData()` (multiple overloads)
 * - `addCellData()` (multiple overloads)
 *
 * \ingroup pythonbindings
 *
 * \tparam Writer The writer class type.
 * \tparam options Additional options for the writer class.
 * \param scope The scope in which to register the class.
 * \param cls The class object to register the methods with.
 */
template <class Writer, class... options>
void registerVtkWriter(pybind11::handle scope, pybind11::class_<Writer, options...> cls) {
  using pybind11::operator""_a;

  using Ikarus::Vtk::DataTag;
  // using Ikarus::Vtk::DataTag::asPointData;

  using Assembler     = typename Writer::Assembler;
  using FE            = typename Writer::FEType;
  using GridView      = typename Writer::GridView;
  using VirtualizedGF = Dune::Vtk::Function<GridView>;

  cls.def(pybind11::init(
              [](std::shared_ptr<Assembler> assembler, Dune::Vtk::FormatTypes format, Dune::Vtk::DataTypes datatype,
                 Dune::Vtk::DataTypes headertype) { return new Writer(assembler, format, datatype, headertype); }),
          pybind11::arg("assembler"), pybind11::arg("format") = Dune::Vtk::FormatTypes::BINARY,
          pybind11::arg("datatype")   = Dune::Vtk::DataTypes::FLOAT32,
          pybind11::arg("headertype") = Dune::Vtk::DataTypes::UINT32, pybind11::keep_alive<1, 2>());

  cls.def("setFormat", &Writer::setFormat);
  cls.def("setDatatype", &Writer::setDatatype);
  cls.def("setHeadertype", &Writer::setHeadertype);

  cls.def("addAllResults", [](Writer& self, DataTag tag) { self.addAllResults(tag); });

  auto addResultImpl = [](Writer& self, const std::string& resType, auto type) {
    bool success = false;
    Dune::Hybrid::forEach(typename FE::SupportedResultTypes(), [&]<typename RT>(RT i) {
      if (resType == toString(i)) {
        success = true;
        self.template addResult<RT::template Rebind>(type);
      }
    });
    if (not success)
      DUNE_THROW(Dune::NotImplemented, "Element " + Dune::className<FE>() + " doesn't support ResultType " + resType);
  };

  cls.def(
      "addResult", [&](Writer& self, const std::string& resType, DataTag tag) { addResultImpl(self, resType, tag); },
      pybind11::arg("resType"), pybind11::arg("tag"));

  cls.def(
      "write", [](Writer& self, const std::string& fileName) { return self.write(fileName); },
      pybind11::arg("fileName"));
}

} // namespace Ikarus::Python

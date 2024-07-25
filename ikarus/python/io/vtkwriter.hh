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
#include <ikarus/io/vtkwriter.hh>

// PYBIND11_MAKE_OPAQUE(std::vector<bool>);
namespace Ikarus::Python {

/**
 * \brief Register Python bindings for a VtkWriter class.
 *
 * \ingroup pythonbindings
 */
template <class Writer, class... options>
void registerVtkWriter(pybind11::handle scope, pybind11::class_<Writer, options...> cls) {
  using pybind11::operator""_a;

  using Ikarus::Vtk::asCellData;
  using Ikarus::Vtk::asPointData;

  using Assembler     = typename Writer::Assembler;
  using FE            = typename Writer::FEType;
  using GridView      = typename Writer::GridView;
  using VirtualizedGF = Dune::Vtk::Function<GridView>;

  cls.def(pybind11::init([](const Assembler& assembler, Dune::Vtk::FormatTypes format, Dune::Vtk::DataTypes datatype,
                            Dune::Vtk::DataTypes headertype) {
            return new Writer(std::make_shared<Assembler>(assembler), format, datatype, headertype);
          }),
          pybind11::arg("assembler"), pybind11::arg("format") = Dune::Vtk::FormatTypes::BINARY,
          pybind11::arg("datatype")   = Dune::Vtk::DataTypes::FLOAT32,
          pybind11::arg("headertype") = Dune::Vtk::DataTypes::UINT32, pybind11::keep_alive<1, 2>());

  cls.def("setFormat", &Writer::setFormat);
  cls.def("setDatatype", &Writer::setDatatype);
  cls.def("setHeadertype", &Writer::setHeadertype);

  cls.def("addAllResultsAsCellData", [](Writer& self) { self.addAllResults(asCellData); });

  cls.def("addAllResultsAsPointData", [](Writer& self) { self.addAllResults(asPointData); });

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
      "addResultAsCellData",
      [&](Writer& self, const std::string& resType) { addResultImpl(self, resType, asCellData); },
      pybind11::arg("resType"));

  cls.def(
      "addResultAsPointData",
      [&](Writer& self, const std::string& resType) { addResultImpl(self, resType, asPointData); },
      pybind11::arg("resType"));

  cls.def(
      "write", [](Writer& self, const std::string& fileName) { return self.write(fileName); },
      pybind11::arg("fileName"));

  // The following bindings are directly copied from /dune/dune-vtk/dune/python/vtk/writer.hh, maybe there is a better
  // way to do this
  cls.def(
      "addPointData",
      [](Writer& writer, VirtualizedGF& f, Dune::Vtk::RangeTypes range, Dune::Vtk::DataTypes data) {
        f.setRangeType(range);
        f.setDataType(data);
        writer.addPointData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("range") = Dune::Vtk::RangeTypes::AUTO,
      pybind11::arg("data") = Dune::Vtk::DataTypes::FLOAT32);
  cls.def(
      "addPointData",
      [](Writer& writer, VirtualizedGF& f, std::string& name, Dune::Vtk::RangeTypes range, Dune::Vtk::DataTypes data) {
        f.setName(name);
        f.setRangeType(range);
        f.setDataType(data);
        writer.addPointData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("name"),
      pybind11::arg("range") = Dune::Vtk::RangeTypes::AUTO, pybind11::arg("data") = Dune::Vtk::DataTypes::FLOAT32);
  cls.def(
      "addPointData",
      [](Writer& writer, VirtualizedGF& f, std::string& name, std::vector<int>& components, Dune::Vtk::RangeTypes range,
         Dune::Vtk::DataTypes data) {
        f.setName(name);
        f.setRangeType(range);
        f.setDataType(data);
        f.setComponents(components);
        writer.addPointData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("name"), pybind11::arg("components"),
      pybind11::arg("range") = Dune::Vtk::RangeTypes::AUTO, pybind11::arg("data") = Dune::Vtk::DataTypes::FLOAT32);
  cls.def(
      "addPointData",
      [](Writer& writer, VirtualizedGF& f, Dune::Vtk::FieldInfo& info) {
        f.setFieldInfo(info);
        writer.addPointData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("info"));
  cls.def(
      "addPointData",
      [](Writer& writer, VirtualizedGF& f, std::vector<int>& components, Dune::Vtk::FieldInfo& info) {
        f.setFieldInfo(info);
        f.setComponents(components);
        writer.addPointData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("components"), pybind11::arg("info"));

  cls.def(
      "addCellData",
      [](Writer& writer, VirtualizedGF& f, Dune::Vtk::RangeTypes range, Dune::Vtk::DataTypes data) {
        f.setRangeType(range);
        f.setDataType(data);
        writer.addCellData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("range") = Dune::Vtk::RangeTypes::AUTO,
      pybind11::arg("data") = Dune::Vtk::DataTypes::FLOAT32);
  cls.def(
      "addCellData",
      [](Writer& writer, VirtualizedGF& f, std::string& name, Dune::Vtk::RangeTypes range, Dune::Vtk::DataTypes data) {
        f.setName(name);
        f.setRangeType(range);
        f.setDataType(data);
        writer.addCellData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("name"),
      pybind11::arg("range") = Dune::Vtk::RangeTypes::AUTO, pybind11::arg("data") = Dune::Vtk::DataTypes::FLOAT32);
  cls.def(
      "addCellData",
      [](Writer& writer, VirtualizedGF& f, std::string& name, std::vector<int>& components, Dune::Vtk::RangeTypes range,
         Dune::Vtk::DataTypes data) {
        f.setName(name);
        f.setRangeType(range);
        f.setDataType(data);
        f.setComponents(components);
        writer.addCellData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("name"), pybind11::arg("components"),
      pybind11::arg("range") = Dune::Vtk::RangeTypes::AUTO, pybind11::arg("data") = Dune::Vtk::DataTypes::FLOAT32);
  cls.def(
      "addCellData",
      [](Writer& writer, VirtualizedGF& f, Dune::Vtk::FieldInfo& info) {
        f.setFieldInfo(info);
        writer.addCellData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("info"));
  cls.def(
      "addCellData",
      [](Writer& writer, VirtualizedGF& f, std::vector<int>& components, Dune::Vtk::FieldInfo& info) {
        f.setFieldInfo(info);
        f.setComponents(components);
        writer.addCellData(f);
      },
      pybind11::keep_alive<1, 2>(), pybind11::arg("f"), pybind11::arg("components"), pybind11::arg("info"));
}

} // namespace Ikarus::Python

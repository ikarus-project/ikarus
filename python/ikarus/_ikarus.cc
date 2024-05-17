// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <dune/common/float_cmp.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/python/finiteelements/material.hh>
#include <ikarus/python/finiteelements/valuewrapper.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
// clang-format off
#define ENUM_BINDINGS(Type)                        \
  py::enum_<Type> Type##Enum(m, #Type);           \
  Type Type##EnumV = Type::BEGIN;                          \
  Ikarus::increment(Type##EnumV);                          \
  for (; Type##EnumV != Type::END; Ikarus::increment(Type##EnumV)) \
    Type##Enum.value(toString(Type##EnumV).c_str(), Type##EnumV);
// clang-format on

PYBIND11_MODULE(_ikarus, m) {
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;
  using namespace Eigen;

  ENUM_BINDINGS(ScalarAffordance);
  ENUM_BINDINGS(VectorAffordance);
  ENUM_BINDINGS(MatrixAffordance);
  ENUM_BINDINGS(FESolutions);
  ENUM_BINDINGS(FEParameter);
  ENUM_BINDINGS(SolverTypeTag);
  ENUM_BINDINGS(EnforcingDBCOption);

  py::class_<std::tuple<ScalarAffordance, VectorAffordance, MatrixAffordance>>(m, "Base");
  auto affordanceCollections3 =
      Dune::Python::insertClass<AffordanceCollection<ScalarAffordance, VectorAffordance, MatrixAffordance>,
                                std::tuple<ScalarAffordance, VectorAffordance, MatrixAffordance>>(
          m, "AffordanceCollection",
          Dune::Python::GenerateTypeName("AffordanceCollection<ScalarAffordance,VectorAffordance,MatrixAffordance>"),
          Dune::Python::IncludeFiles("ikarus/finiteelements/ferequirements.hh"))
          .first;

  affordanceCollections3.def(py::init<ScalarAffordance, VectorAffordance, MatrixAffordance>());
  affordanceCollections3.def_property_readonly_static("elastoStatics",
                                                      [](py::object) { return AffordanceCollections::elastoStatics; });

  using VWd = ValueWrapper<double>;
  auto valueWrapperDouble =
      Dune::Python::insertClass<VWd>(m, "ValueWrapper", Dune::Python::GenerateTypeName("ValueWrapper<double>")).first;
  valueWrapperDouble.def(py::init<double>());
  valueWrapperDouble.def("__repr__", [](const VWd& d) { return std::to_string(d.val); });
  valueWrapperDouble.def("__eq__", [](const VWd& x, const VWd& y) { return Dune::FloatCmp::eq(x.val, y.val); });
  valueWrapperDouble.def("__eq__", [](const VWd& x, const double& y) { return Dune::FloatCmp::eq(x.val, y); });
  valueWrapperDouble.def("__eq__", [](const double& x, const VWd& y) { return Dune::FloatCmp::eq(x, y.val); });
  valueWrapperDouble.def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self += py::self)
      .def(-py::self)
      .def(py::self *= double())
      .def(double() * py::self)
      .def(py::self * double());

  auto materials = m.def_submodule("materials", "This is the submodule for materials in Ikarus");

  pybind11::class_<LinearElasticity> linElastic(materials, "LinearElasticity");
  Ikarus::Python::registerLinearElasticity(materials, linElastic);

  pybind11::class_<StVenantKirchhoff> svk(materials, "StVenantKirchhoff");
  Ikarus::Python::registerStVenantKirchhoff(materials, svk);

  pybind11::class_<NeoHooke> nh(materials, "NeoHooke");
  Ikarus::Python::registerNeoHooke(materials, nh);
}

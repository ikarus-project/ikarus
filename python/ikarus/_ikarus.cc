// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <dune/common/float_cmp.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/python/finiteelements/materials/material.hh>

// since python does not support passing python float by reference to a double&, we have to wrap everything
// see also https://pybind11.readthedocs.io/en/stable/faq.html#limitations-involving-reference-arguments
template <typename T>
struct ValueWrapper
{
  T val;
  ValueWrapper operator+(const ValueWrapper& v) const { return ValueWrapper{val + v.val}; }
  ValueWrapper operator-(const ValueWrapper& v) const { return ValueWrapper{val - v.val}; }
  ValueWrapper operator-() const { return ValueWrapper{-val}; }
  ValueWrapper operator*(T value) const { return ValueWrapper{val * value}; }
  ValueWrapper& operator+=(const ValueWrapper& v) {
    val += v.val;
    return *this;
  }
  ValueWrapper& operator*=(T v) {
    val *= v;
    return *this;
  }

  friend ValueWrapper operator*(T f, const ValueWrapper& v) { return ValueWrapper{f * v.val}; }
};

PYBIND11_MODULE(_ikarus, m) {
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;
  using namespace Eigen;

  py::enum_<ScalarAffordances> enumSA(m, "ScalarAffordances");
  enumSA.value("noAffordance", ScalarAffordances::noAffordance);
  enumSA.value("mechanicalPotentialEnergy", ScalarAffordances::mechanicalPotentialEnergy);
  enumSA.value("microMagneticPotentialEnergy", ScalarAffordances::microMagneticPotentialEnergy);

  py::enum_<VectorAffordances> enumVA(m, "VectorAffordances");
  enumVA.value("noAffordance", VectorAffordances::noAffordance);
  enumVA.value("forces", VectorAffordances::forces);
  enumVA.value("microMagneticForces", VectorAffordances::microMagneticForces);

  py::enum_<MatrixAffordances> enumMA(m, "MatrixAffordances");
  enumMA.value("noAffordance", MatrixAffordances::noAffordance);
  enumMA.value("geometricstiffness", MatrixAffordances::geometricstiffness);
  enumMA.value("mass", MatrixAffordances::mass);
  enumMA.value("materialstiffness", MatrixAffordances::materialstiffness);
  enumMA.value("stiffness", MatrixAffordances::stiffness);
  enumMA.value("microMagneticHessian", MatrixAffordances::microMagneticHessian);
  enumMA.value("stiffnessdiffBucklingVector", MatrixAffordances::stiffnessdiffBucklingVector);

  py::enum_<FESolutions> feSol(m, "FESolutions");
  feSol.value("noSolution", FESolutions::noSolution);
  feSol.value("displacement", FESolutions::displacement);
  feSol.value("velocity", FESolutions::velocity);
  feSol.value("director", FESolutions::director);
  feSol.value("magnetizationAndVectorPotential", FESolutions::magnetizationAndVectorPotential);

  py::enum_<FEParameter> fePar(m, "FEParameter");
  fePar.value("noParameter", FEParameter::noParameter);
  fePar.value("loadfactor", FEParameter::loadfactor);
  fePar.value("time", FEParameter::time);

  py::enum_<ResultType> resType(m, "ResultType");
  resType.value("noType", ResultType::noType);
  resType.value("magnetization", ResultType::magnetization);
  resType.value("gradientNormOfMagnetization", ResultType::gradientNormOfMagnetization);
  resType.value("vectorPotential", ResultType::vectorPotential);
  resType.value("divergenceOfVectorPotential", ResultType::divergenceOfVectorPotential);
  resType.value("BField", ResultType::BField);
  resType.value("HField", ResultType::HField);
  resType.value("cauchyStress", ResultType::cauchyStress);
  resType.value("PK2Stress", ResultType::PK2Stress);
  resType.value("linearStress", ResultType::linearStress);
  resType.value("director", ResultType::director);

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

  using FEreq   = FERequirements<Ref<VectorXd>>;
  auto includes = Dune::Python::IncludeFiles{"ikarus/finiteelements/ferequirements.hh"};
  auto lv =
      Dune::Python::insertClass<FEreq>(
          m, "FERequirements", Dune::Python::GenerateTypeName("FERequirements<Eigen::Ref<Eigen::VectorXd>>"), includes)
          .first;
  lv.def(py::init());
  lv.def("addAffordance", [](FEreq& self, const ScalarAffordances& affordances) { self.addAffordance(affordances); });
  lv.def("addAffordance", [](FEreq& self, const VectorAffordances& affordances) { self.addAffordance(affordances); });
  lv.def("addAffordance", [](FEreq& self, const MatrixAffordances& affordances) { self.addAffordance(affordances); });
  lv.def(
      "insertGlobalSolution",
      [](FEreq& self, FESolutions solType, Ref<VectorXd> solVec) {
        self.insertGlobalSolution(std::move(solType), solVec);
      },
      "solutionType"_a, "solutionVector"_a.noconvert());
  lv.def(
      "getGlobalSolution", [](FEreq& self, FESolutions solType) { return self.getGlobalSolution(std::move(solType)); },
      py::return_value_policy::reference_internal);
  lv.def(
      "insertParameter",
      [](FEreq& self, FEParameter parType, ValueWrapper<double>& parVal) {
        self.insertParameter(std::move(parType), parVal.val);
      },
      py::keep_alive<1, 3>(), "FEParameter"_a, "parameterValue"_a.noconvert());

  lv.def("getParameter", [](const FEreq& self, FEParameter parType) { return self.getParameter(std::move(parType)); });

  auto materials = m.def_submodule("materials", "This is the submodule for materials in Ikarus");

  pybind11::class_<LinearElasticity> linElastic(materials, "LinearElasticity");
  Ikarus::Python::registerLinearElasticity(materials, linElastic);

  pybind11::class_<StVenantKirchhoff> svk(materials, "StVenantKirchhoff");
  Ikarus::Python::registerStVenantKirchhoff(materials, svk);

  pybind11::class_<NeoHooke> nh(materials, "NeoHooke");
  Ikarus::Python::registerNeoHooke(materials, nh);
}

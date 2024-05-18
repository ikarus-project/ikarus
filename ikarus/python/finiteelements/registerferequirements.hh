// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/python/pybind11/pybind11.h>

#include <ikarus/python/finiteelements/valuewrapper.hh>

namespace Ikarus::Python {

template <class FE, class... options>
void registerFERequirement(pybind11::handle scope, pybind11::class_<FE, options...> cls) {
  using FERequirements     = typename FE::Requirement;
  using SolutionVectorType = typename FERequirements::SolutionVectorType;
  using ParameterType      = typename FERequirements::ParameterType;

  using ParameterType = typename FERequirements::ParameterType;
  using namespace pybind11::literals;

  cls.def(
      "createRequirement", [](pybind11::object /* self */) { return FERequirements(); },
      pybind11::return_value_policy::copy);

  auto includes              = Dune::Python::IncludeFiles{"ikarus/finiteelements/ferequirements.hh"};
  auto [lv, isNotRegistered] = Dune::Python::insertClass<FERequirements>(
      scope, "FERequirements", Dune::Python::GenerateTypeName(Dune::className<FERequirements>()), includes);
  if (isNotRegistered) {
    lv.def(pybind11::init());
    lv.def(pybind11::init<SolutionVectorType&, ParameterType&>());

    lv.def(
        "insertGlobalSolution",
        [](FERequirements& self, SolutionVectorType solVec) { self.insertGlobalSolution(solVec); },
        "solutionVector"_a.noconvert());
    lv.def(
        "globalSolution", [](FERequirements& self) { return self.globalSolution(); },
        pybind11::return_value_policy::reference_internal);
    lv.def(
        "insertParameter", [](FERequirements& self, ValueWrapper<double>& parVal) { self.insertParameter(parVal.val); },
        pybind11::keep_alive<1, 2>(), "parameterValue"_a.noconvert());

    lv.def("parameter", [](const FERequirements& self) { return self.parameter(); });
  }
}
} // namespace Ikarus::Python
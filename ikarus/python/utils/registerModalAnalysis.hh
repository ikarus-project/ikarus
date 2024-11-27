// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/bitsetvector.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#define EIGEN_DEFAULT_TO_ROW_MAJOR 1

namespace Ikarus::Python {

template <class ModalAnalysis, class... options>
void registerModalAnalysis(pybind11::handle scope, pybind11::class_<ModalAnalysis, options...> cls) {
  using pybind11::operator""_a;

  using Assembler           = typename ModalAnalysis::Assembler;
  using DirichletValuesType = typename Assembler::DirichletValuesType;
  using FEContainer         = typename ModalAnalysis::FEContainer;

  cls.def(pybind11::init([](const pybind11::list& fes, const DirichletValuesType& dirichletValues) {
            FEContainer fesV = fes.template cast<FEContainer>();
            return new ModalAnalysis(std::move(fesV), dirichletValues);
          }),
          pybind11::keep_alive<1, 3>());

  cls.def("compute", &ModalAnalysis::compute, pybind11::return_value_policy::copy);
  cls.def("angularFrequencies", &ModalAnalysis::angularFrequencies, pybind11::return_value_policy::copy);
  cls.def("naturalFrequencies", &ModalAnalysis::naturalFrequencies, pybind11::return_value_policy::copy);
  cls.def("squaredAngularFrequencies", &ModalAnalysis::squaredAngularFrequencies, pybind11::return_value_policy::copy);
  cls.def("eigenmodes", &ModalAnalysis::eigenmodes, pybind11::return_value_policy::copy);
  cls.def("nev", &ModalAnalysis::nev);
  // cls.def()
}

} // namespace Ikarus::Python

// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/io/resultfunction.hh>
#include <dune/python/pybind11/pybind11.h>


namespace Ikarus::Python {
  /*
   * For now no user function
   */

  template <class ResultFunction, class... options>
  void registerResultFunction(pybind11::handle scope, pybind11::class_<ResultFunction, options...> cls) {
    using pybind11::operator""_a;

    using ElementType = typename ResultFunction::ElementType;
    using ResultRequirements = typename ResultFunction::ResultRequirements;
    using GridView = typename ResultFunction::GridView;
    using ctype = typename ResultFunction::ctype;

    // Was für keep alive, was für nurse was für patient
    cls.def(pybind11::init([](std::vector<ElementType>* fes, const ResultRequirements& req) {
                  return new ResultFunction(fes, req);
                }), pybind11::keep_alive<1, 2>());


  }
}

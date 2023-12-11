// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

  // Python wrapper for the FVAssembler C++ class
  template <class Basis, class... options>
  void registerBasis(pybind11::handle scope, pybind11::class_<Basis, options...> cls) {
    using pybind11::operator""_a;

    using GridView       = typename Basis::GridView;
    using PreBasis       = typename Basis::PreBasis;
    using UntouchedBasis = typename Basis::UntouchedBasis;
    using FlatBasis      = typename Basis::FlatBasis;

    cls.def(pybind11::init([](const UntouchedBasis& gb) { return new Basis(gb.preBasis()); }));

    cls.def(
        "flat", [](Basis& self) { return self.flat(); }, pybind11::return_value_policy::reference);

    cls.def(
        "untouched", [](Basis& self) { return self.untouched(); }, pybind11::return_value_policy::reference);
  }

}  // namespace Ikarus::Python

// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteElements/feRequirements.hh>

namespace Ikarus::Python {

  // Python wrapper for the FVAssembler C++ class
  template <class BoundaryPatch, class... options>
  void registerBoundaryPatch(pybind11::handle scope, pybind11::class_<BoundaryPatch, options...> cls) {
    using pybind11::operator""_a;

    using GridView = typename BoundaryPatch::GridView;

    cls.def(pybind11::init([](const GridView& gv, Eigen::Ref<Eigen::VectorX<bool>> vec) {
              Dune::BitSetVector<1> bitSetVector;
              bitSetVector.resize(vec.size());
              for (size_t i = 0; i < vec.size(); ++i)
                bitSetVector[i] = vec[i];
              return new BoundaryPatch(gv, bitSetVector);
            }),
            pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());
  }

}  // namespace Ikarus::Python

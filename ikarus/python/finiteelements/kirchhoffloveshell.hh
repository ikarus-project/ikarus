// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "linearelastic.hh"

#include <dune/fufem/boundarypatch.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/mechanics/kirchhoffloveshell.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

  template <class KirchhoffLoveShell, class... options>
  void registerKirchhoffLoveShell(pybind11::handle scope, pybind11::class_<KirchhoffLoveShell, options...> cls) {
    registerElement<false, KirchhoffLoveShell, options...>(scope, cls);
    using GlobalBasis    = typename KirchhoffLoveShell::Basis;
    using FlatBasis      = typename KirchhoffLoveShell::FlatBasis;
    using GridView       = typename GlobalBasis::GridView;
    using Element        = typename KirchhoffLoveShell::Element;
    using Traits         = typename KirchhoffLoveShell::Traits;
    using FErequirements = typename KirchhoffLoveShell::FERequirementType;

    using LoadFunction = std::function<Eigen::Vector<double, Traits::worlddim>(Eigen::Vector<double, Traits::worlddim>,
                                                                               const double&)>;
    cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu,
                              double thickness, const LoadFunction volumeLoad) {
              return new KirchhoffLoveShell(basis, element, emod, nu, thickness, volumeLoad);
            }),
            pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

    cls.def(
        pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu, double thickness) {
          return new KirchhoffLoveShell(basis, element, emod, nu, thickness);
        }),
        pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

    cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu,
                              double thickness, const LoadFunction volumeLoad, const BoundaryPatch<GridView>& bp,
                              const LoadFunction neumannBoundaryLoad) {
              return new KirchhoffLoveShell(basis, element, emod, nu, thickness, volumeLoad, &bp, neumannBoundaryLoad);
            }),
            pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>(), pybind11::keep_alive<1, 8>());
  }

}  // namespace Ikarus::Python

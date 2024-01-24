// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file kirchhoffloveshell.hh
 * \brief Python bindings for the Kirchhoff-Love shell element
 */

#pragma once

#include "registerelement.hh"

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

/**
 * \brief Register Python bindings for a KirchhoffLoveShell class.
 *
 * This function registers Python bindings for a KirchhoffLoveShell class, allowing it to be used in Python scripts.
 * The registered class will have several initializers with different sets of parameters.
 *
 * \tparam KirchhoffLoveShell The KirchhoffLoveShell class to be registered.
 * \tparam options Variadic template parameters for additional options when defining the Python class.
 *
 * \param scope A Pybind11 handle representing the Python scope where the class should be registered.
 * \param cls The Pybind11 class template to be used for registering the KirchhoffLoveShell class.
 *
 * \ingroup pythonbindings
 */
template <class KirchhoffLoveShell, class... options>
void registerKirchhoffLoveShell(pybind11::handle scope, pybind11::class_<KirchhoffLoveShell, options...> cls) {
  registerElement<false, KirchhoffLoveShell, options...>(scope, cls);
  using GlobalBasis    = typename KirchhoffLoveShell::Basis;
  using FlatBasis      = typename KirchhoffLoveShell::FlatBasis;
  using GridView       = typename GlobalBasis::GridView;
  using Element        = typename KirchhoffLoveShell::Element;
  using Traits         = typename KirchhoffLoveShell::Traits;
  using FERequirements = typename KirchhoffLoveShell::FERequirementType;

  using LoadFunction = std::function<Eigen::Vector<double, Traits::worlddim>(
      Dune::FieldVector<double, Traits::worlddim>, const double&)>;
  cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu, double thickness,
                            const LoadFunction volumeLoad) {
            return new KirchhoffLoveShell(basis, element, emod, nu, thickness, volumeLoad);
          }),
          pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

  cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu,
                            double thickness) { return new KirchhoffLoveShell(basis, element, emod, nu, thickness); }),
          pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

  cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu, double thickness,
                            const LoadFunction volumeLoad, const BoundaryPatch<GridView>& bp,
                            const LoadFunction neumannBoundaryLoad) {
            return new KirchhoffLoveShell(basis, element, emod, nu, thickness, volumeLoad, &bp, neumannBoundaryLoad);
          }),
          pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>(), pybind11::keep_alive<1, 8>());
}

} // namespace Ikarus::Python

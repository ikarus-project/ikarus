// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file dirichletvalues.hh
 * \brief Python bindings for DirichletValues
 */

#pragma once

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/stl_bind.h>

#include <ikarus/finiteelements/ferequirements.hh>
// PYBIND11_MAKE_OPAQUE(std::vector<bool>);
namespace Ikarus::Python {

/**
 * \brief Register Python bindings for a DirichletValues class.
 *
 * This function registers Python bindings for a DirichletValues class, allowing it to be used in Python scripts.
 * The registered class will have an initializer that takes a `Basis` object. It exposes several member functions to
 * Python:
 *   - `fixBoundaryDOFs(f)`: Fixes boundary degrees of freedom using a user-defined function `f`.
 *   - `fixBoundaryDOFsUsingLocalView(f)`: Fixes boundary degrees of freedom using a user-defined function `f` with a
 * `LocalView` argument.
 *   - `fixBoundaryDOFsUsingLocalViewAndIntersection(f)`: Fixes boundary degrees of freedom using a user-defined
 * function `f` with `LocalView` and `Intersection` arguments.
 *   - `fixDOFs(f)`: Fixes boundary degrees of freedom using a user-defined function `f` with the boolean vector and
 * the basis as arguments.
 *
 * \tparam DirichletValues The DirichletValues class to be registered.
 * \tparam options Variadic template parameters for additional options when defining the Python class.
 *
 * \param scope A Pybind11 handle representing the Python scope where the class should be registered.
 * \param cls The Pybind11 class template to be used for registering the DirichletValues class.
 *
 * \ingroup pythonbindings
 */
template <class DirichletValues, class... options>
void registerDirichletValues(pybind11::handle scope, pybind11::class_<DirichletValues, options...> cls) {
  using pybind11::operator""_a;

  using Basis        = typename DirichletValues::Basis;
  using BackendType  = typename DirichletValues::BackendType;
  using FlagsType    = typename DirichletValues::FlagsType;
  using MultiIndex   = typename Basis::MultiIndex;
  using LocalView    = typename Basis::LocalView;
  using Intersection = typename Basis::GridView::Intersection;

  pybind11::module scopedf = pybind11::module::import("dune.functions");
  typedef Dune::Python::LocalViewWrapper<Basis> LocalViewWrapper;
  auto includes = Dune::Python::IncludeFiles{"dune/python/functions/globalbasis.hh"};
  auto lv       = Dune::Python::insertClass<LocalViewWrapper>(
                scopedf, "LocalView",
                Dune::Python::GenerateTypeName("Dune::Python::LocalViewWrapper", Dune::MetaType<Basis>()), includes)
                .first;

  cls.def(pybind11::init([](const Basis& basis) { return new DirichletValues(basis); }), pybind11::keep_alive<1, 2>());

  // Eigen::Ref needed due to https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html#pass-by-reference
  cls.def("fixBoundaryDOFs",
          [](DirichletValues& self, const std::function<void(Eigen::Ref<Eigen::VectorX<bool>>, int)>& f) {
            auto lambda = [&](BackendType& vec, const MultiIndex& indexGlobal) {
              // we explicitly only allow flat indices
              f(vec.vector(), indexGlobal[0]);
            };
            self.fixBoundaryDOFs(lambda);
          });

  cls.def("fixBoundaryDOFsUsingLocalView",
          [](DirichletValues& self,
             const std::function<void(Eigen::Ref<Eigen::VectorX<bool>>, int, LocalViewWrapper&)>& f) {
            auto lambda = [&](BackendType& vec, int localIndex, LocalView& lv) {
              auto lvWrapper = LocalViewWrapper(lv.globalBasis());
              // this can be simplified when
              // https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/418 becomes available
              pybind11::object obj = pybind11::cast(lv.element());
              lvWrapper.bind(obj);
              f(vec.vector(), localIndex, lvWrapper);
            };
            self.fixBoundaryDOFs(lambda);
          });

  cls.def(
      "fixBoundaryDOFsUsingLocalViewAndIntersection",
      [](DirichletValues& self,
         const std::function<void(Eigen::Ref<Eigen::VectorX<bool>>, int, LocalViewWrapper&, const Intersection&)>& f) {
        auto lambda = [&](BackendType& vec, int localIndex, LocalView& lv, const Intersection& intersection) {
          auto lvWrapper = LocalViewWrapper(lv.globalBasis());
          // this can be simplified when
          // https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/418 becomes available
          pybind11::object obj = pybind11::cast(lv.element());
          lvWrapper.bind(obj);
          f(vec.vector(), localIndex, lvWrapper, intersection);
        };
        self.fixBoundaryDOFs(lambda);
      });

  cls.def("fixDOFs",
          [](DirichletValues& self, const std::function<void(const Basis&, Eigen::Ref<Eigen::VectorX<bool>>)>& f) {
            auto lambda = [&](const Basis& basis, BackendType& vec) {
              // we explicitly only allow flat indices
              f(basis, vec.vector());
            };
            self.fixDOFs(lambda);
          });
}

} // namespace Ikarus::Python

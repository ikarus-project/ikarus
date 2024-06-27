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

template <class DirichletValues, typename CppVisitor>
void forwardCorrectFunction(DirichletValues& dirichletValues, const pybind11::function& functor,
                            CppVisitor&& cppFunction) {
  using Basis        = typename DirichletValues::Basis;
  using Intersection = typename Basis::GridView::Intersection;
  using BackendType  = typename DirichletValues::BackendType;
  using MultiIndex   = typename Basis::MultiIndex;

  using LocalViewWrapper = Dune::Python::LocalViewWrapper<Basis>;

  using FixBoundaryDOFsWithGlobalIndexFunction = std::function<void(Eigen::Ref<Eigen::VectorX<bool>>, int)>;
  using FixBoundaryDOFsWithLocalViewFunction =
      std::function<void(Eigen::Ref<Eigen::VectorX<bool>>, int, LocalViewWrapper&)>;
  using FixBoundaryDOFsWithIntersectionFunction =
      std::function<void(Eigen::Ref<Eigen::VectorX<bool>>, int, LocalViewWrapper&, const Intersection&)>;

  // Disambiguate by number of arguments, as casting doesn't properly work with functions
  pybind11::module inspect_module = pybind11::module::import("inspect");
  pybind11::object result         = inspect_module.attr("signature")(functor).attr("parameters");
  size_t numParams                = pybind11::len(result);

  if (numParams == 2) {
    auto& function = functor.template cast<const FixBoundaryDOFsWithGlobalIndexFunction>();
    auto lambda    = [&](BackendType& vec, const MultiIndex& indexGlobal) {
      // we explicitly only allow flat indices
      function(vec.vector(), indexGlobal[0]);
    };
    cppFunction(lambda);

  } else if (numParams == 3) {
    auto& function = functor.template cast<const FixBoundaryDOFsWithLocalViewFunction>();
    auto lambda    = [&](BackendType& vec, int localIndex, auto& lv) {
      auto lvWrapper = LocalViewWrapper(lv.rootLocalView());
      function(vec.vector(), localIndex, lvWrapper);
    };
    cppFunction(lambda);

  } else if (numParams == 4) {
    auto& function = functor.template cast<const FixBoundaryDOFsWithIntersectionFunction>();
    auto lambda    = [&](BackendType& vec, int localIndex, auto& lv, const Intersection& intersection) {
      auto lvWrapper = LocalViewWrapper(lv.rootLocalView());
      function(vec.vector(), localIndex, lvWrapper, intersection);
    };
    cppFunction(lambda);

  } else {
    DUNE_THROW(Dune::NotImplemented, "fixBoundaryDOFs: A function with this signature is not supported");
  }
}

/**
 * \brief Register Python bindings for a DirichletValues class.
 *
 * This function registers Python bindings for a DirichletValues class, allowing it to be used in Python scripts.
 * The registered class will have an initializer that takes a `Basis` object. It exposes several member functions to
 * Python:
 *   - `fixBoundaryDOFs(f)`: Fixes boundary degrees of freedom. This function can be overloaded with the following
 * arguments:
 *    - using a user-defined function `f`.
 *    -  using a user-defined function `f` with a `LocalView` argument.
 *    -  using a user-defined function `f` with `LocalView` and `Intersection` arguments.
 *    - `fixDOFs(f)`: Fixes boundary degrees of freedom using a user-defined function `f` with the boolean vector and
 * the basis as arguments.
 *  -  - `fixBoundaryDOFsOfSubSpaceBasis(f)`: Same implementation but you can pass a index of a child basis of a Power
 * Basis
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
  auto lv_      = Dune::Python::insertClass<LocalViewWrapper>(
                 scopedf, "LocalView",
                 Dune::Python::GenerateTypeName("Dune::Python::LocalViewWrapper", Dune::MetaType<Basis>()), includes)
                 .first;

  cls.def(pybind11::init([](const Basis& basis) { return new DirichletValues(basis); }), pybind11::keep_alive<1, 2>());

  cls.def_property_readonly("container", &DirichletValues::container);
  cls.def_property_readonly("size", &DirichletValues::size);
  cls.def_property_readonly("fixedDOFsize", &DirichletValues::fixedDOFsize);
  cls.def("isConstrained", [](DirichletValues& self, std::size_t i) -> bool { return self.isConstrained(i); });

  cls.def("fixIthDOF", [](DirichletValues& self, std::size_t i) { self.fixIthDOF(MultiIndex(std::array{i})); });

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

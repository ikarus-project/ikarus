// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file registerElement.hh
 * \brief Python bindings for a generic finite element
 */

#pragma once

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
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

/**
 * \brief Register Python bindings for a generic finite element class.
 *
 * This function registers Python bindings for a generic finite element class, allowing it to be used in Python
 * scripts. The registered class will have multiple initializers with different sets of parameters and member
 * functions to calculate results and access properties.
 *
 * \tparam defaultInitializers If true, include default initializers for the finite element class.
 * \tparam FE The generic finite element class to be registered.
 * \tparam options Variadic template parameters for additional options when defining the Python class.
 *
 * \param scope A Pybind11 handle representing the Python scope where the class should be registered.
 * \param cls The Pybind11 class template to be used for registering the finite element class.
 *
 * \ingroup pythonbindings
 */
template <bool defaultInitializers = true, class FE, class... options>
void registerElement(pybind11::handle scope, pybind11::class_<FE, options...> cls) {
  using pybind11::operator""_a;

  using GlobalBasis    = typename FE::Basis;
  using FlatBasis      = typename FE::FlatBasis;
  using GridView       = typename GlobalBasis::GridView;
  using Element        = typename FE::Element;
  using Traits         = typename FE::Traits;
  using FERequirements = typename FE::FERequirementType;

  if constexpr (defaultInitializers)
    cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu) {
              return new FE(basis, element, emod, nu);
            }),
            pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

  using LoadFunction = std::function<Eigen::Vector<double, Traits::worlddim>(
      Dune::FieldVector<double, Traits::worlddim>, const double&)>;
  if constexpr (defaultInitializers)
    cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu,
                              const LoadFunction volumeLoad) { return new FE(basis, element, emod, nu, volumeLoad); }),
            pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());
  if constexpr (defaultInitializers)
    cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu,
                              const LoadFunction volumeLoad, const BoundaryPatch<GridView>& bp,
                              const LoadFunction neumannBoundaryLoad) {
              return new FE(basis, element, emod, nu, volumeLoad, &bp, neumannBoundaryLoad);
            }),
            pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>(), pybind11::keep_alive<1, 7>());

  pybind11::module scopedf = pybind11::module::import("dune.functions");

  typedef Dune::Python::LocalViewWrapper<FlatBasis> LocalViewWrapper;
  auto includes = Dune::Python::IncludeFiles{"dune/python/functions/globalbasis.hh"};
  auto lv       = Dune::Python::insertClass<LocalViewWrapper>(
                scopedf, "LocalViewWrapper",
                Dune::Python::GenerateTypeName("Dune::Python::LocalViewWrapperWrapper", Dune::MetaType<FlatBasis>()),
                includes)
                .first;
  lv.def("bind", &LocalViewWrapper::bind);
  lv.def("unbind", &LocalViewWrapper::unbind);
  lv.def("index", [](const LocalViewWrapper& localView, int index) { return localView.index(index); });
  lv.def("__len__", [](LocalViewWrapper& self) -> int { return self.size(); });

  Dune::Python::Functions::registerTree<typename LocalViewWrapper::Tree>(lv);
  lv.def("tree", [](const LocalViewWrapper& view) { return view.tree(); });

  auto basisName = Dune::className<FlatBasis>();
  auto entry     = Dune::Python::insertClass<FlatBasis>(scopedf, "DefaultGlobalBasis", pybind11::buffer_protocol(),
                                                    Dune::Python::GenerateTypeName(basisName),
                                                    Dune::Python::IncludeFiles{"dune/python/functions/globalbasis.hh"});
  if (entry.second)
    Dune::Python::registerGlobalBasis(scopedf, entry.first);

  cls.def(
      "localView",
      [](FE& self) -> LocalViewWrapper {
        auto lvWrapped = LocalViewWrapper(self.localView().globalBasis());
        // this can be simplified when https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/418
        // becomes available
        pybind11::object obj = pybind11::cast(self.localView().element());
        lvWrapped.bind(obj);
        return lvWrapped;
      },
      pybind11::keep_alive<0, 1>());
  cls.def("calculateScalar", [](FE& self, const FERequirements& req) { return self.calculateScalar(req); });
  cls.def("calculateVector", [](FE& self, const FERequirements& req, Eigen::Ref<Eigen::VectorXd> vec) {
    return self.calculateVector(req, vec);
  });
  cls.def(
      "calculateMatrix",
      [](FE& self, const FERequirements& req, Eigen::Ref<Eigen::MatrixXd> mat) {
        return self.calculateMatrix(req, mat);
      },
      pybind11::arg("FERequirements"), pybind11::arg("elementMatrix").noconvert());

  if constexpr (requires { std::declval<FE>().materialTangent(); })
    cls.def("materialTangent", [](FE& self) { return self.materialTangent(); });

  cls.def(
      "calculateAt",
      [](FE& self, const FERequirements& req, const Dune::FieldVector<double, Traits::mydim>& local,
         ResultType resType) {
        if (resType == ResultType::linearStress)
          return self.template calculateAt<ResultType::linearStress>(req, local);
        else
          DUNE_THROW(Dune::NotImplemented, "Linear-lastic element only supports linearStress as result.");
      },
      pybind11::arg("feRequirements"), pybind11::arg("local"), pybind11::arg("resultType"));
}

} // namespace Ikarus::Python

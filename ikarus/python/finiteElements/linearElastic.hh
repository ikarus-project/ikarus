// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

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

#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

  template <class LinearElastic, class... options>
  void registerLinearElastic(pybind11::handle scope, pybind11::class_<LinearElastic, options...> cls) {
    using pybind11::operator""_a;

    using GlobalBasis    = typename LinearElastic::Basis;
    using FlatBasis      = typename LinearElastic::FlatBasis;
    using GridView       = typename GlobalBasis::GridView;
    using Element        = typename LinearElastic::Element;
    using Traits         = typename LinearElastic::Traits;
    using FErequirements = typename LinearElastic::FERequirementType;

    cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu) {
              return new LinearElastic(basis, element, emod, nu);
            }),
            pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

    using LoadFunction = std::function<Eigen::Vector<double, Traits::worlddim>(Eigen::Vector<double, Traits::worlddim>,
                                                                               const double&)>;
    cls.def(pybind11::init(
                [](const GlobalBasis& basis, const Element& element, double emod, double nu,
                   const LoadFunction volumeLoad) { return new LinearElastic(basis, element, emod, nu, volumeLoad); }),
            pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 3>());

    cls.def(pybind11::init([](const GlobalBasis& basis, const Element& element, double emod, double nu,
                              const LoadFunction volumeLoad, const BoundaryPatch<GridView>& bp,
                              const LoadFunction neumannBoundaryLoad) {
              return new LinearElastic(basis, element, emod, nu, volumeLoad, &bp, neumannBoundaryLoad);
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
    auto entry     = Dune::Python::insertClass<FlatBasis>(
        scopedf, "DefaultGlobalBasis", pybind11::buffer_protocol(), Dune::Python::GenerateTypeName(basisName),
        Dune::Python::IncludeFiles{"dune/python/functions/globalbasis.hh"});
    if (entry.second) Dune::Python::registerGlobalBasis(scopedf, entry.first);

    cls.def(
        "localView",
        [](LinearElastic& self) -> LocalViewWrapper {
          auto lvWrapped = LocalViewWrapper(self.localView().globalBasis());
          // this can be simplified when https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/418
          // becomes available
          pybind11::object obj = pybind11::cast(self.localView().element());
          lvWrapped.bind(obj);
          return lvWrapped;
        },
        pybind11::keep_alive<0, 1>());
    cls.def("calculateScalar",
            [](LinearElastic& self, const FErequirements& req) { return self.calculateScalar(req); });
    cls.def("calculateVector", [](LinearElastic& self, const FErequirements& req, Eigen::Ref<Eigen::VectorXd> vec) {
      return self.calculateVector(req, vec);
    });
    cls.def(
        "calculateMatrix",
        [](LinearElastic& self, const FErequirements& req, Eigen::Ref<Eigen::MatrixXd> mat) {
          return self.calculateMatrix(req, mat);
        },
        pybind11::arg("FErequirements"), pybind11::arg("elementMatrix").noconvert());

    cls.def("getMaterialTangent", [](LinearElastic& self) { return self.getMaterialTangent(); });
  }

}  // namespace Ikarus::Python

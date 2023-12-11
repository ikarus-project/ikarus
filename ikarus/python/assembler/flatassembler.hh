// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/stl_bind.h>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

#define MAKE_ASSEMBLER_REGISTERY_FUNCTION(name)                                                                                \
  template <class Assembler, class... options>                                                                                 \
  void register##name(pybind11::handle scope, pybind11::class_<Assembler, options...> cls) {                                   \
    using pybind11::operator""_a;                                                                                              \
    using FEContainer         = typename Assembler::FEContainer;                                                               \
    using Basis               = typename Assembler::Basis;                                                                     \
    using DirichletValuesType = typename Assembler::DirichletValuesType;                                                       \
    using FERequirementType   = typename Assembler::FERequirementType;                                                         \
    pybind11::module m        = pybind11::module::import("ikarus");                                                            \
    cls.def(pybind11::init([](const pybind11::list& fes, const DirichletValuesType& dirichletValues) {                         \
              /*here a copy of the whole vector of fes takes place! There is no way to prevent this if we want that            \
               * the user can pass native python lists here, see                                                               \
               * https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html */                                           \
              FEContainer fesV = fes.template cast<FEContainer>();                                                             \
              return new Assembler(std::move(fesV), dirichletValues);                                                          \
            }),                                                                                                                \
            pybind11::keep_alive<1, 3>());                                                                                     \
                                                                                                                               \
    /* sparse matrices need to be copied to python therefore we remove the reference of the return type, see */                \
    /*  https://github.com/pybind/pybind11/blob/cbb876cc7b02c5f57e715cbc2c46ead3d1fbcf79/tests/test_eigen_matrix.cpp#L332-L341 \
     */                                                                                                                        \
    cls.def(                                                                                                                   \
        "getMatrix",                                                                                                           \
        [](Assembler& self, const FERequirementType& req) -> std::remove_cvref_t<decltype(self.getMatrix(req))> {              \
          return self.getMatrix(req);                                                                                          \
        },                                                                                                                     \
        pybind11::return_value_policy::copy);                                                                                  \
    cls.def(                                                                                                                   \
        "getReducedMatrix",                                                                                                    \
        [](Assembler& self, const FERequirementType& req)                                                                      \
            -> std::remove_cvref_t<decltype(self.getReducedMatrix(req))> { return self.getReducedMatrix(req); },               \
        pybind11::return_value_policy::copy);                                                                                  \
    cls.def(                                                                                                                   \
        "getVector", [](Assembler& self, const FERequirementType& req) { return self.getVector(req); },                        \
        pybind11::return_value_policy::reference);                                                                             \
    cls.def(                                                                                                                   \
        "getScalar", [](Assembler& self, const FERequirementType& req) { return self.getScalar(req); },                        \
        pybind11::return_value_policy::copy);                                                                                  \
    cls.def(                                                                                                                   \
        "getReducedVector", [](Assembler& self, const FERequirementType& req) { return self.getReducedVector(req); },          \
        pybind11::return_value_policy::reference);                                                                             \
    cls.def(                                                                                                                   \
        "createFullVector",                                                                                                    \
        [](Assembler& self, Eigen::Ref<const Eigen::VectorXd> redVec) { return self.createFullVector(redVec); },               \
        pybind11::return_value_policy::move);                                                                                  \
    cls.def(                                                                                                                   \
        "reducedSize", [](Assembler& self) { return self.reducedSize(); }, pybind11::return_value_policy::copy);               \
  }

  MAKE_ASSEMBLER_REGISTERY_FUNCTION(SparseFlatAssembler);
  MAKE_ASSEMBLER_REGISTERY_FUNCTION(DenseFlatAssembler);

}  // namespace Ikarus::Python

// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file flatassembler.hh
 * \brief Python bindings for assemblers
 */

#pragma once

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/stl_bind.h>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

#define MAKE_ASSEMBLER_REGISTERY_FUNCTION(name)                                                                                \
  /**                                                                                                                          \
   \brief Register Python bindings for a name class.                                                                           \
\n                                                                                                                             \
    \details This function registers Python bindings for a name class, allowing it to be used in Python scripts.   \n          \
  This function is a result of the macro `MAKE_ASSEMBLER_REGISTERY_FUNCTION(name)`.   \n                                       \
  \n                                                                                                                           \
   The registered class will have an initializer that takes a list of finite elements (`fes`) and a                            \
`DirichletValuesType` object. \n It exposes several member functions to Python:  \n <ul> <li> `matrix(req)`:                   \
Returns a dense matrix based on the specified `FERequirementType`.    \n <li> `vector(req,affo,qt)`: Returns a vector          \
based on the specified `FERequirementType`.   \n <li> `scalar(req,affo)`: Returns a scalar based on the specified              \
`FERequirementType`.    \n <li> `createFullVector(redVec)`: Creates a full vector from a reduced vector.   \n <li>             \
`reducedSize()`: Returns the size of the reduced space.   \n                                                                   \
  </ul> \n                                                                                                                     \
   \tparam Assembler The name class to be registered.                                                                          \
   \tparam options Variadic template parameters for additional options when defining the Python class.                         \
   \param scope A Pybind11 handle representing the Python scope where the class should be registered.                          \
   \param cls The Pybind11 class template to be used for registering the name class.                                           \
   \ingroup pythonbindings                                                                                                     \
  */                                                                                                                           \
  template <class Assembler, class... options>                                                                                 \
  void register##name(pybind11::handle scope, pybind11::class_<Assembler, options...> cls) {                                   \
    using pybind11::operator""_a;                                                                                              \
    using FEContainer              = typename Assembler::FEContainer;                                                          \
    using Basis                    = typename Assembler::Basis;                                                                \
    using DirichletValuesType      = typename Assembler::DirichletValuesType;                                                  \
    using AffordanceCollectionType = typename Assembler::AffordanceCollectionType;                                             \
    using FERequirementType        = typename Assembler::FERequirement;                                                        \
    pybind11::module m             = pybind11::module::import("ikarus");                                                       \
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
        "matrix",                                                                                                              \
        [](Assembler& self, const FERequirementType& req, Ikarus::MatrixAffordance affo,                                       \
           Ikarus::EnforcingDBCOption qt) -> std::remove_cvref_t<decltype(self.getMatrix(req, affo))> {                        \
          return self.matrix(req, affo, qt);                                                                                   \
        },                                                                                                                     \
        pybind11::return_value_policy::copy);                                                                                  \
                                                                                                                               \
    cls.def(                                                                                                                   \
        "matrix", [](Assembler& self) -> std::remove_cvref_t<decltype(self.getMatrix())> { return self.matrix(); },            \
        pybind11::return_value_policy::copy);                                                                                  \
                                                                                                                               \
    cls.def(                                                                                                                   \
        "matrix",                                                                                                              \
        [](Assembler& self, Ikarus::EnforcingDBCOption qt) -> std::remove_cvref_t<decltype(self.getMatrix(qt))> {              \
          return self.matrix(qt);                                                                                              \
        },                                                                                                                     \
        pybind11::return_value_policy::copy);                                                                                  \
                                                                                                                               \
    cls.def(                                                                                                                   \
        "vector",                                                                                                              \
        [](Assembler& self, const FERequirementType& req, Ikarus::VectorAffordance affo,                                       \
           Ikarus::EnforcingDBCOption qt) { return self.vector(req, affo, qt); },                                              \
        pybind11::return_value_policy::reference);                                                                             \
                                                                                                                               \
    cls.def("vector", [](Assembler& self) { return self.vector(); }, pybind11::return_value_policy::reference);                \
                                                                                                                               \
    cls.def(                                                                                                                   \
        "vector", [](Assembler& self, Ikarus::EnforcingDBCOption qt) { return self.vector(qt); },                              \
        pybind11::return_value_policy::reference);                                                                             \
                                                                                                                               \
    cls.def(                                                                                                                   \
        "scalar",                                                                                                              \
        [](Assembler& self, const FERequirementType& req, Ikarus::ScalarAffordance affo) {                                     \
          return self.scalar(req, affo);                                                                                       \
        },                                                                                                                     \
        pybind11::return_value_policy::copy);                                                                                  \
                                                                                                                               \
    cls.def("scalar", [](Assembler& self) { return self.scalar(); }, pybind11::return_value_policy::copy);                     \
                                                                                                                               \
    cls.def(                                                                                                                   \
        "createFullVector",                                                                                                    \
        [](Assembler& self, Eigen::Ref<const Eigen::VectorXd> redVec) { return self.createFullVector(redVec); },               \
        pybind11::return_value_policy::move);                                                                                  \
    cls.def("reducedSize", [](Assembler& self) { return self.reducedSize(); }, pybind11::return_value_policy::copy);           \
    cls.def("bind", [](Assembler& self, const FERequirementType& req, AffordanceCollectionType affo,                           \
                       EnforcingDBCOption qt = EnforcingDBCOption::Full) { return self.bind(req, affo, qt); });                \
    cls.def("bind", [](Assembler& self, const FERequirementType& req) { return self.bind(req); });                             \
    cls.def("bind", [](Assembler& self, const AffordanceCollectionType affo) { return self.bind(affo); });                     \
    cls.def("bind", [](Assembler& self, const EnforcingDBCOption qt) { return self.bind(qt); });                               \
    cls.def("bound", &Assembler::bound);                                                                                       \
    cls.def("boundToRequirement", &Assembler::boundToRequirement);                                                             \
    cls.def("boundToAffordanceCollection", &Assembler::boundToAffordanceCollection);                                           \
    cls.def("boundToEnforcingDBCOption", &Assembler::boundToEnforcingDBCOption);                                               \
    cls.def("requirement", &Assembler::requirement);                                                                           \
    cls.def("affordanceCollection", &Assembler::affordanceCollection);                                                         \
    cls.def("enforcingDBCOption", &Assembler::enforcingDBCOption);                                                             \
  }
// TODO Alex add binding functions
MAKE_ASSEMBLER_REGISTERY_FUNCTION(SparseFlatAssembler);
MAKE_ASSEMBLER_REGISTERY_FUNCTION(DenseFlatAssembler);

} // namespace Ikarus::Python

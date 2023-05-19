// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/classname.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteElements/mechanics/materials.hh>

#define MAKE_MaterialFunction(clsName, materialName, functionname, vecSize)                                            \
  clsName.def(                                                                                                         \
      #functionname,                                                                                                   \
      [](materialName& self, const std::string& straintag, Eigen::Ref<const Eigen::Vector<double, vecSize>> eVoigt_) { \
        Eigen::Vector<double, vecSize> eVoigt = eVoigt_;                                                               \
        if (straintag == toString(StrainTags::rightCauchyGreenTensor))                                                 \
          return self.template functionname<StrainTags::rightCauchyGreenTensor>(eVoigt);                               \
        else if (straintag == toString(StrainTags::greenLagrangian))                                                   \
          return self.template functionname<StrainTags::greenLagrangian>(eVoigt);                                      \
        else if (straintag == toString(StrainTags::linear))                                                            \
          DUNE_THROW(Dune::MathError, "Passing linear strain to " + std::string(#materialName)                         \
                                          + " does not makes sense use LinearElastic class");                          \
        else if (straintag == toString(StrainTags::displacementGradient))                                              \
          DUNE_THROW(Dune::MathError,                                                                                  \
                     "Passing displacementGradient strain in 6d Voigt notation does not make any sense!");             \
        else if (straintag == toString(StrainTags::deformationGradient))                                               \
          DUNE_THROW(Dune::MathError,                                                                                  \
                     "Passing deformationGradient strain in 6d Voigt notation does not make any sense!");              \
      },                                                                                                               \
      "StrainName"_a, "strainVector"_a);

#define MAKE_MATERIAL_REGISTERY_FUNCTION(Materialname, vecSize)                                                                \
  template <class Materialname, class... options>                                                                              \
  void register##Materialname(pybind11::handle scope, pybind11::class_<Materialname, options...> cls##Materialname) {          \
    using pybind11::operator""_a;                                                                                              \
    namespace py = pybind11;                                                                                                   \
    cls##Materialname.def(pybind11::init([](double emod, double nu) {                                                          \
                            auto matParameter                                                                                  \
                                = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = emod, .nu = nu});                    \
                            return new Materialname(matParameter);                                                             \
                          }),                                                                                                  \
                          "Young's modulus"_a, "Poisson's ratio"_a);                                                           \
    MAKE_MaterialFunction(cls##Materialname, Materialname, storedEnergy, vecSize);                                             \
    MAKE_MaterialFunction(cls##Materialname, Materialname, stresses, vecSize);                                                 \
    MAKE_MaterialFunction(cls##Materialname, Materialname, tangentModuli, vecSize);                                            \
                                                                                                                               \
    using PlaneStressClass = decltype(planeStress(std::declval<Materialname>()));                                              \
    auto includes          = Dune::Python::IncludeFiles{"ikarus/finiteElements/mechanics/materials.hh"};                       \
    auto pS                = Dune::Python::insertClass<PlaneStressClass>(                                                      \
                  scope, "PlaneStressClass",                                                                    \
                  Dune::Python::GenerateTypeName(                                                               \
                                     "Ikarus::VanishingStress<std::array<Ikarus::Impl::StressIndexPair, "                      \
                                                    "3ul>{{Ikarus::Impl::StressIndexPair{2ul, 1ul}, Ikarus::Impl::StressIndexPair{2ul,0ul}, " \
                                                    "Ikarus::Impl::StressIndexPair{2ul, 2ul}}},"                                              \
                                     + Dune::className<Materialname>() + ">"),                                                 \
                  includes)                                                                                     \
                  .first;                                                                                                      \
    MAKE_MaterialFunction(pS, PlaneStressClass, storedEnergy, 3);                                                              \
    MAKE_MaterialFunction(pS, PlaneStressClass, stresses, 3);                                                                  \
    MAKE_MaterialFunction(pS, PlaneStressClass, tangentModuli, 3);                                                             \
    MAKE_MaterialFunction(pS, PlaneStressClass, storedEnergy, 6);                                                              \
    MAKE_MaterialFunction(pS, PlaneStressClass, stresses, 6);                                                                  \
    MAKE_MaterialFunction(pS, PlaneStressClass, tangentModuli, 6);                                                             \
    cls##Materialname.def("asPlaneStress", [](Materialname& self) {                                                            \
      return planeStress(self);                                                                                                \
    }); /* no keep_alive since planeStress copies the material */                                                              \
    cls##Materialname.def("asShellMaterial", [](Materialname& self) {                                                          \
      return shellMaterial(self);                                                                                              \
    }); /* no keep_alive since shellMaterial copies the material */                                                            \
    cls##Materialname.def("asbeamMaterial", [](Materialname& self) {                                                           \
      return beamMaterial(self);                                                                                               \
    }); /* no keep_alive since beamMaterial copies the material */                                                             \
  }

#define MAKE_MATERIAL_CLASS_IN_MODULE(Materialname, args)                                                      \
  auto includes##Materialname = Dune::Python::IncludeFiles{"ikarus/finiteElements/mechanics/materials.hh"};    \
  auto cls##Materialname                                                                                       \
      = Dune::Python::insertClass<Ikarus::Materialname<args>>(                                                 \
            m, #Materialname, Dune::Python::GenerateTypeName("Ikarus::" + std::string(#Materialname<##args>)), \
            includes##Materialname)                                                                            \
            .first;                                                                                            \
  cls##Materialname.def(pybind11::init([](double emod, double nu) {                                            \
                          auto matParameter                                                                    \
                              = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = emod, .nu = nu});      \
                          return new Materialname(matParameter);                                               \
                        }),                                                                                    \
                        "Young's modulus"_a, "Poisson's ratio"_a);                                             \
  MAKE_MaterialFunction(Materialname<##args>, storedEnergy);                                                   \
  MAKE_MaterialFunction(Materialname<##args>, stresses);                                                       \
  MAKE_MaterialFunction(Materialname<##args>, tangentModuli);

namespace Ikarus::Python {

  MAKE_MATERIAL_REGISTERY_FUNCTION(StVenantKirchhoff, 6);
  MAKE_MATERIAL_REGISTERY_FUNCTION(NeoHooke, 6);
}  // namespace Ikarus::Python

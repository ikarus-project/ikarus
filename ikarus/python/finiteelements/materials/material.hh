// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/classname.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/utils/concepts.hh>

#define MAKE_MaterialFunction(clsName, materialName, functionname, vecSize)                                            \
  clsName.def(                                                                                                         \
      #functionname,                                                                                                   \
      [](materialName& self, const std::string& straintag, Eigen::Ref<const Eigen::Vector<double, vecSize>> eVoigt_) { \
        if constexpr (not Concepts::IsMaterial<LinearElasticityT, materialName>) {                                     \
          Eigen::Vector<double, vecSize> eVoigt = eVoigt_;                                                             \
          if (straintag == toString(StrainTags::rightCauchyGreenTensor))                                               \
            return self.template functionname<StrainTags::rightCauchyGreenTensor>(eVoigt);                             \
          else if (straintag == toString(StrainTags::greenLagrangian))                                                 \
            return self.template functionname<StrainTags::greenLagrangian>(eVoigt);                                    \
          else if (straintag == toString(StrainTags::linear))                                                          \
            DUNE_THROW(Dune::MathError, "Passing linear strain to " + std::string(#materialName)                       \
                                            + " does not makes sense use LinearElastic class");                        \
          else if (straintag == toString(StrainTags::displacementGradient))                                            \
            DUNE_THROW(Dune::MathError,                                                                                \
                       "Passing displacementGradient strain in 6d Voigt notation does not make any sense!");           \
          else if (straintag == toString(StrainTags::deformationGradient))                                             \
            DUNE_THROW(Dune::MathError,                                                                                \
                       "Passing deformationGradient strain in 6d Voigt notation does not make any sense!");            \
          else                                                                                                         \
            DUNE_THROW(Dune::MathError, straintag + "is not a valid strain tag.");                                     \
        } else {                                                                                                       \
          Eigen::Vector<double, vecSize> eVoigt = eVoigt_; /* Linear elastic path */                                   \
          if (straintag == toString(StrainTags::linear))                                                               \
            return self.template functionname<StrainTags::linear>(eVoigt);                                             \
          else                                                                                                         \
            DUNE_THROW(Dune::MathError, "Linear elastic material only accepts linear strains!");                       \
        }                                                                                                              \
        __builtin_unreachable();                                                                                       \
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
                          "Material constructor that takes Young's modulus E and Poisson's ratio nu", "E"_a, "nu"_a);          \
    MAKE_MaterialFunction(cls##Materialname, Materialname, storedEnergy, vecSize);                                             \
    MAKE_MaterialFunction(cls##Materialname, Materialname, stresses, vecSize);                                                 \
    MAKE_MaterialFunction(cls##Materialname, Materialname, tangentModuli, vecSize);                                            \
                                                                                                                               \
    using PlaneStressClass = decltype(planeStress(std::declval<Materialname>()));                                              \
    auto includes          = Dune::Python::IncludeFiles{"ikarus/finiteelements/mechanics/materials.hh"};                       \
    auto pS                = Dune::Python::insertClass<PlaneStressClass>(                                                      \
                  scope, std::string("PlaneStress_") + #Materialname,                                           \
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
    using shellMaterialClass = decltype(shellMaterial(std::declval<Materialname>()));                                          \
    auto shellmaterial                                                                                                         \
        = Dune::Python::insertClass<shellMaterialClass>(                                                                       \
              scope, std::string("Shell_") + #Materialname,                                                                    \
              Dune::Python::GenerateTypeName("Ikarus::VanishingStress<std::array<Ikarus::Impl::StressIndexPair, "              \
                                             "1ul>{{Ikarus::Impl::StressIndexPair{2ul, 2ul}}},"                                \
                                             + Dune::className<Materialname>() + ">"),                                         \
              includes)                                                                                                        \
              .first;                                                                                                          \
    MAKE_MaterialFunction(shellmaterial, shellMaterialClass, storedEnergy, 5);                                                 \
    MAKE_MaterialFunction(shellmaterial, shellMaterialClass, stresses, 5);                                                     \
    MAKE_MaterialFunction(shellmaterial, shellMaterialClass, tangentModuli, 5);                                                \
    MAKE_MaterialFunction(shellmaterial, shellMaterialClass, storedEnergy, 6);                                                 \
    MAKE_MaterialFunction(shellmaterial, shellMaterialClass, stresses, 6);                                                     \
    MAKE_MaterialFunction(shellmaterial, shellMaterialClass, tangentModuli, 6);                                                \
    cls##Materialname.def("asShellMaterial", [](Materialname& self) {                                                          \
      return shellMaterial(self);                                                                                              \
    }); /* no keep_alive since shellMaterial copies the material */                                                            \
    using beamMaterialClass = decltype(beamMaterial(std::declval<Materialname>()));                                            \
    auto beammaterial       = Dune::Python::insertClass<beamMaterialClass>(                                                    \
                            scope, std::string("Beam_") + #Materialname,                                                 \
                            Dune::Python::GenerateTypeName(                                                              \
                                      "Ikarus::VanishingStress<std::array<Ikarus::Impl::StressIndexPair, "                     \
                                            "2ul>{{Impl::StressIndexPair{1, 1},Ikarus::Impl::StressIndexPair{2ul, 2ul}}},"           \
                                      + Dune::className<Materialname>() + ">"),                                                \
                            includes)                                                                                    \
                            .first;                                                                                            \
    MAKE_MaterialFunction(beammaterial, beamMaterialClass, storedEnergy, 4);                                                   \
    MAKE_MaterialFunction(beammaterial, beamMaterialClass, stresses, 4);                                                       \
    MAKE_MaterialFunction(beammaterial, beamMaterialClass, tangentModuli, 4);                                                  \
    MAKE_MaterialFunction(beammaterial, beamMaterialClass, storedEnergy, 6);                                                   \
    MAKE_MaterialFunction(beammaterial, beamMaterialClass, stresses, 6);                                                       \
    MAKE_MaterialFunction(beammaterial, beamMaterialClass, tangentModuli, 6);                                                  \
    cls##Materialname.def("asBeamMaterial", [](Materialname& self) {                                                           \
      return beamMaterial(self);                                                                                               \
    }); /* no keep_alive since beamMaterial copies the material */                                                             \
  }

#define MAKE_MATERIAL_CLASS_IN_MODULE(Materialname, args)                                                      \
  auto includes##Materialname = Dune::Python::IncludeFiles{"ikarus/finiteelements/mechanics/materials.hh"};    \
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
                        "Material constructor that takes Young's modulus E and Poisson's ratio nu"             \
                        "E"_a,                                                                                 \
                        "nu"_a);                                                                               \
  MAKE_MaterialFunction(Materialname<##args>, storedEnergy);                                                   \
  MAKE_MaterialFunction(Materialname<##args>, stresses);                                                       \
  MAKE_MaterialFunction(Materialname<##args>, tangentModuli);

namespace Ikarus::Python {

  MAKE_MATERIAL_REGISTERY_FUNCTION(LinearElasticity, 6);
  MAKE_MATERIAL_REGISTERY_FUNCTION(StVenantKirchhoff, 6);
  MAKE_MATERIAL_REGISTERY_FUNCTION(NeoHooke, 6);
}  // namespace Ikarus::Python

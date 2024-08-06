// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file material.hh
 * \brief Python bindings for materials
 */

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
            DUNE_THROW(Dune::MathError, "Passing linear strain to " + std::string(#materialName) +                     \
                                            " does not makes sense use LinearElastic class");                          \
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

namespace Ikarus::Python {

template <class Material, size_t vecSize, class... options>
void registerMaterial(pybind11::handle scope, pybind11::class_<Material, options...> cls) {
  using pybind11::operator""_a;

  std::string materialname = Material::name();
  cls.def(pybind11::init([](double emod, double nu) {
            auto matParameter = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = emod, .nu = nu});
            return new Material(matParameter);
          }),
          "Material constructor that takes Young's modulus E and Poisson's ratio nu", "E"_a, "nu"_a);
  MAKE_MaterialFunction(cls, Material, storedEnergy, vecSize);
  MAKE_MaterialFunction(cls, Material, stresses, vecSize);
  MAKE_MaterialFunction(cls, Material, tangentModuli, vecSize);

  using PlaneStressClass = decltype(planeStress(std::declval<Material>()));
  auto includes          = Dune::Python::IncludeFiles{"ikarus/finiteelements/mechanics/materials.hh"};
  auto pS                = Dune::Python::insertClass<PlaneStressClass>(
                scope, std::string("PlaneStress_") + materialname,
                Dune::Python::GenerateTypeName(
                    "Ikarus::VanishingStress<std::array<Ikarus::Impl::StressIndexPair, "
                                   "3ul>{{Ikarus::Impl::StressIndexPair{2ul, 1ul}, Ikarus::Impl::StressIndexPair{2ul,0ul}, "
                                   "Ikarus::Impl::StressIndexPair{2ul, 2ul}}}," +
                    Dune::className<Material>() + ">"),
                includes)
                .first;
  MAKE_MaterialFunction(pS, PlaneStressClass, storedEnergy, 3);
  MAKE_MaterialFunction(pS, PlaneStressClass, stresses, 3);
  MAKE_MaterialFunction(pS, PlaneStressClass, tangentModuli, 3);
  MAKE_MaterialFunction(pS, PlaneStressClass, storedEnergy, 6);
  MAKE_MaterialFunction(pS, PlaneStressClass, stresses, 6);
  MAKE_MaterialFunction(pS, PlaneStressClass, tangentModuli, 6);
  cls.def("asPlaneStress",
          [](Material& self) { return planeStress(self); }); /* no keep_alive since planeStress copies the material */

  using PlaneStrainClass = decltype(planeStrain(std::declval<Material>()));
  auto pStrain =
      Dune::Python::insertClass<PlaneStrainClass>(
          scope, std::string("PlaneStrain_") + materialname,
          Dune::Python::GenerateTypeName("Ikarus::PlaneStrain<" + Dune::className<Material>() + ">"), includes)
          .first;
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, storedEnergy, 3);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, stresses, 3);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, tangentModuli, 3);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, storedEnergy, 6);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, stresses, 6);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, tangentModuli, 6);

  cls.def("asPlaneStrain",
          [](Material& self) { return planeStrain(self); }); /* no keep_alive since planeStress copies the material */
  using ShellMaterialClass = decltype(shellMaterial(std::declval<Material>()));
  auto shellmaterial =
      Dune::Python::insertClass<ShellMaterialClass>(
          scope, std::string("Shell_") + materialname,
          Dune::Python::GenerateTypeName("Ikarus::VanishingStress<std::array<Ikarus::Impl::StressIndexPair, "
                                         "1ul>{{Ikarus::Impl::StressIndexPair{2ul, 2ul}}}," +
                                         Dune::className<Material>() + ">"),
          includes)
          .first;
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, storedEnergy, 5);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, stresses, 5);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, tangentModuli, 5);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, storedEnergy, 6);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, stresses, 6);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, tangentModuli, 6);

  cls.def("asShellMaterial", [](Material& self) {
    return shellMaterial(self);
  }); /* no keep_alive since shellMaterial copies the material */
  using BeamMaterialClass = decltype(beamMaterial(std::declval<Material>()));
  auto beammaterial       = Dune::Python::insertClass<BeamMaterialClass>(
                          scope, std::string("Beam_") + materialname,
                          Dune::Python::GenerateTypeName(
                              "Ikarus::VanishingStress<std::array<Ikarus::Impl::StressIndexPair, "
                                    "2ul>{{Impl::StressIndexPair{1, 1},Ikarus::Impl::StressIndexPair{2ul, 2ul}}}," +
                              Dune::className<Material>() + ">"),
                          includes)
                          .first;
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, storedEnergy, 4);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, stresses, 4);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, tangentModuli, 4);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, storedEnergy, 6);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, stresses, 6);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, tangentModuli, 6);
  cls.def("asBeamMaterial",
          [](Material& self) { return beamMaterial(self); }); /* no keep_alive since beamMaterial copies the material */
}

#define MAKE_MATERIAL_REGISTERY_FUNCTION(name, vecSize)                                     \
  template <class Material, class... options>                                               \
  void register##name(pybind11::handle scope, pybind11::class_<Material, options...> cls) { \
    Ikarus::Python::registerMaterial<Material, vecSize>(scope, cls);                        \
  }

MAKE_MATERIAL_REGISTERY_FUNCTION(LinearElasticity, 6);
MAKE_MATERIAL_REGISTERY_FUNCTION(StVenantKirchhoff, 6);
MAKE_MATERIAL_REGISTERY_FUNCTION(NeoHooke, 6);
} // namespace Ikarus::Python

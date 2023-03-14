// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteElements/mechanics/materials/strainConversions.hh>
#include <ikarus/finiteElements/mechanics/materials/tags.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {
  template <class MaterialImpl>
  struct Material;

  template <auto stressIndexPair, typename MaterialImpl>
  struct VanishingStress;

  template <class MaterialImpl_>
  struct Material {
    using MaterialImpl = MaterialImpl_;

    static constexpr bool isReduced
        = Std::IsSpecializationNonTypeAndTypes<Ikarus::VanishingStress, MaterialImpl>::value;

    /* Const accessor to the underlying material   */
    constexpr MaterialImpl const &impl() const  // CRTP
    {
      return static_cast<MaterialImpl const &>(*this);
    }

    /* Const accessor to the underlying material   */
    constexpr MaterialImpl &impl()  // CRTP
    {
      return static_cast<MaterialImpl &>(*this);
    }

    /* Name of the material    */
    [[nodiscard]] std::string name() const { return impl().nameImpl(); }

    /**
     * Return the stored potential energy of the material
     * @tparam The strain tag, which indicates, which strain tensor components are passed
     * @tparam Derived The underlying Eigen type
     * @param Eraw The strain tensor components, which can be passed in Voigt notation or matrix notation
     * @return Scalar return of stored energy
     */
    template <StrainTags tag, typename Derived>
    [[nodiscard]] auto storedEnergy(const Eigen::MatrixBase<Derived> &Eraw) const {
      decltype(auto) Ev = enlargeIfReduced<Material>(Eraw);
      decltype(auto) E  = transformStrain<tag, MaterialImpl::strainTag>(Ev);
      if constexpr (Concepts::EigenVector<Derived>) {  // receiving vector means voigt notation
        if constexpr (MaterialImpl::energyAcceptsVoigt)
          return impl().storedEnergyImpl(toVoigt(E));
        else
          return impl().storedEnergyImpl(E);
      } else
        return impl().storedEnergyImpl(E);
    }

    /**
     * The stresses of the material
     * @tparam The strain tag, which indicates, which strain tensor components are passed
     * @tparam Derived The underlying Eigen type
     * @param Eraw The strain tensor components, which can be passed in Voigt notation or matrix notation
     * @return Vectorial or Matrix return of stresses
     */
    template <StrainTags tag, bool voigt = true, typename Derived>
    [[nodiscard]] auto stresses(const Eigen::MatrixBase<Derived> &Eraw) const {
      decltype(auto) Ev = enlargeIfReduced<Material>(Eraw);
      decltype(auto) E  = transformStrain<tag, MaterialImpl::strainTag>(Ev);
      if constexpr (voigt and MaterialImpl::stressToVoigt == false)
        // user request a Voigt shaped return but material is not able to. Therefore, we transform it here.
        return toVoigt(stressesMaybeTransformInputToVoigt<false>(E), false);
      else
        return stressesMaybeTransformInputToVoigt<voigt>(E);
    }

    /**
     * The tangentModuli of the material
     * @tparam The strain tag, which indicates, which strain tensor components are passed
     * @tparam Derived The underlying Eigen type
     * @param Eraw The strain tensor components, which can be passed in Voigt notation or matrix notation
     * @return tangent moduli in voigt notation or as fourth order tensor
     */
    template <StrainTags tag, bool voigt = true, typename Derived>
    [[nodiscard]] auto tangentModuli(const Eigen::MatrixBase<Derived> &Eraw) const {
      decltype(auto) Ev = enlargeIfReduced<Material>(Eraw);
      decltype(auto) E  = transformStrain<tag, MaterialImpl::strainTag>(Ev);
      if constexpr (voigt and MaterialImpl::moduliToVoigt == false)
        // user request a Voigt shaped return but material is not able to. Therefore, we transform it here.
        return toVoigt(tangentModuliMaybeTransformInputToVoigt<false>(E));
      else
        return tangentModuliMaybeTransformInputToVoigt<voigt>(E);
    }

    /* Rebind material to different scalar type (Useful for automatic differentiation  */
    template <typename ScalarTypeOther>
    auto rebind() const {
      return impl().template rebind<ScalarTypeOther>();
    }

  private:
    template <bool voigt = true, typename Derived>
    auto stressesMaybeTransformInputToVoigt(const Eigen::MatrixBase<Derived> &E) const {
      if constexpr (Concepts::EigenVector<Derived>) {  // receiving vector means Voigt notation
        if constexpr (MaterialImpl::stressAcceptsVoigt)
          return impl().template stressesImpl<voigt>(E);
        else  // material is not able to accept Voigt shaped Input. Therefore, we transform it before.
          return impl().template stressesImpl<voigt>(fromVoigt(E.derived()));
      } else
        return impl().template stressesImpl<voigt>(E.derived());
    }

    template <bool voigt = true, typename Derived>
    auto tangentModuliMaybeTransformInputToVoigt(const Eigen::MatrixBase<Derived> &E) const {
      if constexpr (Concepts::EigenVector<Derived>) {  // receiving vector means voigt notation
        if constexpr (MaterialImpl::moduliAcceptsVoigt)
          return impl().template tangentModuliImpl<voigt>(E);
        else
          return impl().template tangentModuliImpl<voigt>(fromVoigt(E.derived()));
      } else
        return impl().template tangentModuliImpl<voigt>(E.derived());
    }
  };

}  // namespace Ikarus

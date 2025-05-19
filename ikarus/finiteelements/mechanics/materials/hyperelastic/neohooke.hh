// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file neohooke.hh
 * \brief Implementation of the NeoHooke material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/utils/lambertw.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Neo-Hookean material model.
* \ingroup materials
*  The energy is computed as
*  \f[ \psi(\BC) = \frac{\mu}{2} (\tr \BC-3- 2 \log \sqrt{\det \BC}) + \frac{\la}{2} (\log \sqrt{\det \BC})^2 ,\f]
* where \f$ \BC \f$ denotes the right Cauchy-Green strain tensor.
*
*  The second Piola-Kirchhoff stresses are computed as
*      \f[ \BS(\BC) =\fracpt{\psi(\BC)}{\BC} = \mu (\BI-\BC^{-1}) + \la \log \sqrt{\det \BC}  \BC^{-1},\f]
*
* and the material tangent moduli are computed as
*      \f[ \BBC(\BC) =\fracpt{^2\psi(\BC)}{\BC^2} =  \la \BC^{-1} \otimes  \BC^{-1} + 2 (\mu- \la \log
\sqrt{\det \BC} ) \CI,\f]
*      where \f$ \CI_{IJKL} =  \frac{1}{2}({(\BC^{-1})}^{IK}{(\BC^{-1})}^{JL}+{(\BC^{-1})}^{IL} {(\BC^{-1})}^{JK}).\f$
*
*  \remark See \cite bonet2008nonlinear, Section 6.4.3 for a discussion of this material
 *\tparam ST The underlying scalar type.
 */
template <typename ST>
struct NeoHookeT : public Material<NeoHookeT<ST>>
{
  using ScalarType         = ST;
  static constexpr int dim = 3;
  using StrainMatrix       = Eigen::Matrix<ScalarType, dim, dim>;
  using StressMatrix       = StrainMatrix;
  using MaterialTensor     = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>>;

  using MaterialParameters = LamesFirstParameterAndShearModulus;

  static constexpr auto strainTag              = StrainTags::rightCauchyGreenTensor;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 2;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "NeoHooke"; }

  /**
   * \brief Constructor for NeoHookeT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit NeoHookeT(const MaterialParameters& mpt)
      : materialParameter_{mpt} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!Concepts::EigenVector<Derived>) {
      const auto traceC = C.trace();
      const auto detC   = C.determinant();
      Impl::checkPositiveOrAbort(detC);
      const auto logdetF = log(sqrt(detC));
      return materialParameter_.mu / 2.0 * (traceC - 3 - 2 * logdetF) +
             materialParameter_.lambda / 2.0 * logdetF * logdetF;
    } else
      static_assert(!Concepts::EigenVector<Derived>,
                    "NeoHooke energy can only be called with a matrix and not a vector in Voigt notation");
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  StressMatrix stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        const auto detC = C.determinant();
        Impl::checkPositiveOrAbort(detC);
        const auto logdetF = log(sqrt(detC));
        const auto invC    = C.inverse().eval();
        return (materialParameter_.mu * (StrainMatrix::Identity() - invC) + materialParameter_.lambda * logdetF * invC)
            .eval();
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "NeoHooke can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "NeoHooke does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return MaterialTensor The tangent moduli.
   */
  template <bool voigt, typename Derived>
  MaterialTensor tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      const auto invC = C.inverse().eval();
      const auto detC = C.determinant();
      Impl::checkPositiveOrAbort(detC);
      const auto logdetF = log(sqrt(detC));
      const auto CTinv   = tensorView(invC, std::array<Eigen::Index, 2>({3, 3}));

      MaterialTensor moduli = (materialParameter_.lambda * dyadic(CTinv, CTinv) +
                               2 * (materialParameter_.mu - materialParameter_.lambda * logdetF) *
                                   symTwoSlots(fourthOrderIKJL(invC, invC), {2, 3}))
                                  .eval();
      return moduli;
    } else
      static_assert(voigt == false, "NeoHooke does not support returning tangent moduli in Voigt notation");
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return NeoHookeT<ScalarTypeOther> The rebound NeoHooke material.
   */
  template <typename STO>
  auto rebind() const {
    return NeoHookeT<STO>(materialParameter_);
  }

  /**
   * \brief Computes the strain measure and inverse material tangent for a given stress state.
   * \details Analytical solution is taken from R. Pfefferkorn, S. Bieber, B. Oesterle, M. Bischoff, and P. Betsch,
   * “Improving efficiency and robustness of enhanced assumed strain elements for nonlinear problems,” Numerical Meth
   * Engineering, vol. 122, no. 8, pp. 1911–1939, Apr. 2021, doi: 10.1002/nme.6605.
   *
   * \tparam Derived The derived type of the input matrix.
   * \param Sraw the input stress matrix.
   * \return pair of inverse material tangent and right Cauchy-Green strain tensor in voigt notation.
   */
  template <typename Derived>
  auto materialInversionImpl(const Eigen::MatrixBase<Derived>& Sraw) const {
    decltype(auto) S = Impl::maybeFromVoigt(Sraw, false);
    auto C           = Derived::Zero().eval();
    const auto I     = Derived::Identity();

    auto mu     = materialParameter_.mu;
    auto Lambda = materialParameter_.lambda;
    if (Dune::FloatCmp::eq(Lambda, 0.0))
      C = (I - (1.0 / mu) * S).inverse();
    else {
      const auto A = 1.0 / Lambda * (mu * I - S);
      const auto a = A.determinant();
      checkInvertabilityOrAbort(a);

      const auto b   = 2.0 / dim * mu / Lambda - log(static_cast<double>(dim) / 2) + 1.0 / dim * log(abs(a));
      const auto lnJ = mu / Lambda - static_cast<double>(dim) / 2 * util::lambertW0(Dune::sign(a) * exp(b));
      C              = (mu / Lambda - lnJ) * A.inverse();
    }
    const auto tangentModulus = toVoigt(tangentModuliImpl<false>(C));
    const auto D              = tangentModulus.inverse().eval();
    return std::make_pair(D, toVoigt(C));
  }

private:
  MaterialParameters materialParameter_;

  template <typename ScalarType>
  void checkInvertabilityOrAbort(ScalarType a, ScalarType eps = 1e-10) const {
    const bool req1 = Dune::FloatCmp::ne(a, 0.0, eps);
    const bool req2 = [&]() {
      if constexpr (dim == 3)
        return Dune::FloatCmp::gt(
            a,
            -pow(static_cast<double>(dim) / 2, dim) * exp(-dim - 2 * materialParameter_.mu / materialParameter_.lambda),
            eps);
      else if (dim == 2)
        return Dune::FloatCmp::gt(a, 0.0, eps);
    }();
    if (!(req1 && req2)) {
      std::cerr << "NeoHooke: Requirements for material law inversion not met.";
      abort();
    }
  }
};

/**
 * \brief Alias for NeoHookeT with double as the default scalar type.
 */
using NeoHooke = NeoHookeT<double>;

} // namespace Ikarus::Materials

// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file numericalmaterialinversion.hh
 * \brief Implementation of the numerical scheme for material inverison for generic hyperelastic material models.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Computes numerically an approximation of the complementary strain energy function and returns the material
 * tangent and strain tensor for a given stress state.
 * \details  N. Viebahn, J. Schröder, und P. Wriggers, „An extension of assumed stress finite elements to a general
 * hyperelastic framework“, Adv. Model. and Simul. in Eng. Sci., Bd. 6, Nr. 1, S. 9, Dez. 2019,
 * doi: 10.1186/s40323-019-0133-z.
 *
 * \tparam Material the type of the material
 * \tparam Derived the type of the stress matrix
 * \param mat the material
 * \param S stress matrix
 * \param tol tolerance for the Newton-Raphson solver.
 * \param maxIter maximum number of iterations for the Newton-Raphson solver.
 * \return pair of inverse material tangent and the strain tensor in voigt notation.
 */
template <Concepts::Material Material, typename Derived>
auto numericalMaterialInversion(const Material& mat, const Eigen::MatrixBase<Derived>& S,
                                const Eigen::MatrixBase<Derived>& Estart = Derived::Zero().eval(),
                                const double tol = 1e-12, const int maxIter = 20) {
  static_assert(Concepts::EigenMatrix33<decltype(S)> or Concepts::EigenMatrix22<decltype(S)>);
  constexpr int dim      = Derived::CompileTimeTraits::RowsAtCompileTime;
  constexpr int dimVoigt = (dim * (dim + 1)) / 2;

  using ST       = Derived::Scalar;
  using VoigtVec = Eigen::Vector<ST, dimVoigt>;

  VoigtVec Es     = toVoigt(Estart.derived()); // starting value
  VoigtVec Svoigt = toVoigt(S.derived(), false);

  auto r = VoigtVec::Zero().eval();
  auto D = Eigen::Matrix<ST, dimVoigt, dimVoigt>::Zero().eval();

  for (auto i : Dune::range(maxIter)) {
    r = (Svoigt - mat.template stresses<StrainTags::greenLagrangian, true>(Es)).eval();

    if (r.norm() < tol and i > 0)
      break;
    else if (i == maxIter)
      DUNE_THROW(
          Dune::MathError,
          "Numerical material inversion failed to converge within the maximum number of iterations for the material: " +
              mat.name());

    D = mat.template tangentModuli<StrainTags::greenLagrangian, true>(Es).inverse().eval(); // voigt
    Es += (D * r).eval();
  }
  return std::make_pair(D, transformStrain<StrainTags::greenLagrangian, Material::strainTag>(Es).eval());
}

} // namespace Ikarus::Materials

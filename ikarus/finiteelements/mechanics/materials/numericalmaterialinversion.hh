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
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/nonlinearoperator.hh>
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
                                const Eigen::MatrixBase<Derived>& Estart = Derived::Zero().eval(), double tol = 1e-10,
                                int maxIter = 20) {
  static_assert(Concepts::EigenMatrix33<decltype(S)>);
  Derived Es = Estart; // starting value

  auto r = Eigen::Matrix3d::Zero().eval();
  auto D = Eigen::Matrix<double, 6, 6>::Zero().eval();

  for (auto i : Dune::range(maxIter)) {
    r           = (S - mat.template stresses<StrainTags::greenLagrangian, false>(Es)).eval();
    D           = toVoigt(mat.template tangentModuli<StrainTags::greenLagrangian, false>(Es)).inverse().eval(); // voigt
    auto rVoigt = toVoigt(r, false);
    Es += fromVoigt((D * rVoigt).eval());

    if (rVoigt.norm() < tol)
      break;
    else if (i == maxIter)
      DUNE_THROW(Dune::MathError,
                 "Numerical material inversion failed to converge within  the maximum number of iterations");
  }
  return std::make_pair(D, toVoigt(transformStrain<StrainTags::greenLagrangian, Material::strainTag>(Es).eval()));
}

} // namespace Ikarus::Materials

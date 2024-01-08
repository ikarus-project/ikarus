// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <ranges>

#include <dune/common/fvector.hh>
#include <dune/common/overloadset.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>

#include <Eigen/Core>
namespace Ikarus {

  struct DefaultMembraneStrain {
    template <typename Geometry>
    void pre(const Geometry &geo, const auto &uFunction) const {}

    template <typename Geometry>
    auto value(const Dune::FieldVector<double, 2> &gpPos, const Geometry &geo, const auto &uFunction) const
        -> Eigen::Vector3<typename std::remove_cvref_t<decltype(uFunction)>::ctype> {
      using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
      Eigen::Vector3<ScalarType> epsV;
      const auto J = Dune::toEigen(geo.jacobianTransposed(gpPos));
      using namespace Dune;
      using namespace Dune::DerivativeDirections;
      const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
          uFunction.evaluateDerivative(gpPos,  // Here the gpIndex could be passed
                                       Dune::wrt(spatialAll), Dune::on(Dune::DerivativeDirections::referenceElement)));
      const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();

      epsV << J.row(0).dot(gradu.col(0)) + 0.5 * gradu.col(0).squaredNorm(),
          J.row(1).dot(gradu.col(1)) + 0.5 * gradu.col(1).squaredNorm(), j.row(0).dot(j.row(1));
      return epsV;
    }
    template <typename Geometry, typename ScalarType>
    auto derivative(const Dune::FieldVector<double, 2> &gpPos, const Eigen::Matrix<ScalarType, 2, 3> &jcur,
                    const auto &dNAtGp, const Geometry &geo, const auto &uFunction, const auto &localBasis,
                    const int node) const {
      Eigen::Matrix<ScalarType, 3, 3> bop;
      bop.row(0) = jcur.row(0) * dNAtGp(node, 0);
      bop.row(1) = jcur.row(1) * dNAtGp(node, 1);
      bop.row(2) = jcur.row(0) * dNAtGp(node, 1) + jcur.row(1) * dNAtGp(node, 0);

      return bop;
    }
    template <typename Geometry, typename ScalarType>
    auto secondDerivative(const Dune::FieldVector<double, 2> &gpPos, const auto &dNAtGp, const Geometry &geo,
                          const auto &uFunction, const auto &localBasis, const Eigen::Vector3<ScalarType> &S, int I,
                          int J) const {
      const auto &dN1i                   = dNAtGp(I, 0);
      const auto &dN1j                   = dNAtGp(J, 0);
      const auto &dN2i                   = dNAtGp(I, 1);
      const auto &dN2j                   = dNAtGp(J, 1);
      const ScalarType NS                = dN1i * dN1j * S[0] + dN2i * dN2j * S[1] + (dN1i * dN2j + dN2i * dN1j) * S[2];
      Eigen::Matrix<ScalarType, 3, 3> kg = Eigen::Matrix<double, 3, 3>::Identity() * NS;
      return kg;
    }
  };

}  // namespace Ikarus

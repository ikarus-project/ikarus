// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file displacementgradienttransposed.hh
 * \brief Implementation of the EAS function, where the transpose of the displacement gradient (pre-multiplied with
 * deformation gradient at the center of the element) is enhanced.
 *
 * \ingroup eas
 */

#pragma once

#include <dune/common/fvector.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>
#include <dune/localfefunctions/expressions/greenLagrangeStrains.hh>

#include <Eigen/Core>

#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS {
/**
 * \brief A struct computing the value, first and second derivatives of the Green-Lagrange strain tensor (in Voigt
 * notation), where the transpose of the displacement gradient (pre-multiplied with deformation gradient at the center
 * of the element) is enhanced.
 *
 * \details See \cite glaserFormulationEnhancedStrain1997a for details.
 */
struct DisplacementGradientTransposed
{
  /**
   * \brief Compute the strain vector at a given integration point or its index.
   *
   * \param geo The geometry object providing the transposed Jacobian.
   * \param uFunction The function representing the displacement field.
   * \param gpPos The position of the integration point.
   * \param easFunction The EAS function.
   * \param alpha The coefficients of the EAS function.
   *
   * \tparam GEO The type of the geometry object.
   * \tparam EAST The type of the EAS function.
   *
   * \return The Green-Lagrange strain vector at the given integration point.
   */
  template <typename GEO, typename EAST>
  static auto value(const GEO& geo, const auto& uFunction, const Dune::FieldVector<double, GEO::mydimension>& gpPos,
                    const EAST& easFunction, const auto& alpha) {
    using ST                                = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    constexpr int myDim                     = GEO::mydimension;
    const Eigen::Matrix<ST, myDim, myDim> H = computeDisplacementGradient(geo, uFunction, gpPos, easFunction, alpha);
    const Eigen::Matrix<ST, myDim, myDim> E = 0.5 * (H + H.transpose() + H.transpose() * H);
    return toVoigt(E, true).eval();
  }

  /**
   * \brief Compute the first derivative of the Green-Lagrange strain w.r.t d or alpha for a given node and integration
   * point.
   *
   * \param geo The geometry object of the finite element.
   * \param uFunction The function representing the displacement field.
   * \param localBasis The local basis object.
   * \param gpPos The position of the integration point.
   * \param gpIndex The index of the integration point.
   * \param node The FE node index (defaults to sNaN).
   * \param alpha The coefficients of the EAS function.
   * \param easFunction The EAS function.
   *
   * \tparam wrtCoeff An integer indicating the coefficient w.r.t which the first derivative of Green-Lagrange strain is
   * to be returned (0 = E,d; 1 = E,a).
   * \tparam GEO The type of the geometry object.
   * \tparam EAST The type of the EAS function.
   *
   * \return The first derivative of the Green-Lagrange strain w.r.t d or alpha.
   */
  template <int wrtCoeff, typename GEO, typename EAST>
  static auto firstDerivative(const GEO& geo, const auto& uFunction, const auto& localBasis, const auto& gpIndex,
                              const Dune::FieldVector<double, GEO::mydimension>& gpPos, const EAST& easFunction,
                              const auto& alpha, const int node = sNaN) {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    using ST                         = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    constexpr int myDim              = GEO::mydimension;
    constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
    static_assert(myDim == 2 or myDim == 3,
                  "An enhancement of displacement gradient (transposed) can only be performed for the 2D or 3D case.");
    constexpr int strainSize = myDim * (myDim + 1) / 2;
    using MatrixType         = Eigen::Matrix<ST, myDim, myDim>;
    const auto centerPos     = centerPosition(geo);

    const MatrixType F = computeDeformationGradient(geo, uFunction, gpPos, easFunction, alpha)
                             .transpose()
                             .eval(); // transposed such that the rows are u_{,1} and u_{,2}
    const typename EAST::HType Harray = easFunction(gpPos);
    const auto g0                     = F.row(0);
    const auto g1                     = F.row(1);

    if constexpr (wrtCoeff == 0) { // E,d
      const auto gradUdI  = uFunction.evaluateDerivative(gpPos, wrt(spatialAll, coeff(node)), on(gridElement));
      const auto gradUdI0 = uFunction.evaluateDerivative(centerPos, wrt(spatialAll, coeff(node)), on(gridElement));
      Eigen::Vector<ST, myDim> dN, dN0, dNtilde;
      dNtilde.setZero();
      for (const auto i : Dune::range(myDim)) {
        dN[i]  = getDiagonalEntry(gradUdI[i], 0);
        dN0[i] = getDiagonalEntry(gradUdI0[i], 0);
      }

      for (const auto i : Dune::range(myDim))
        for (const auto t : Dune::range(enhancedStrainSize))
          for (const auto p : Dune::range(myDim))
            dNtilde[i] += Harray[t](i, p) * alpha[t] * dN0[p];

      dN += dNtilde;
      Eigen::Matrix<ST, strainSize, myDim> bopI;
      bopI.row(0) = dN[0] * g0;
      bopI.row(1) = dN[1] * g1;
      if constexpr (myDim == 2)
        bopI.row(2) = dN[1] * g0 + dN[0] * g1;
      else {
        const auto g2 = F.row(2);
        bopI.row(2)   = dN[2] * g2;
        bopI.row(3)   = dN[2] * g1 + dN[1] * g2;
        bopI.row(4)   = dN[2] * g0 + dN[0] * g2;
        bopI.row(5)   = dN[1] * g0 + dN[0] * g1;
      }
      return bopI.eval();
    } else if constexpr (wrtCoeff == 1) { // E,a
      const MatrixType Fc0 = MatrixType::Identity() + compatibleDisplacementGradient<GEO>(uFunction, centerPos);
      Eigen::Matrix<ST, strainSize, enhancedStrainSize> M;
      M.setZero(); // zeros returned if Harray is empty
      for (const auto q : Dune::range(enhancedStrainSize)) {
        const MatrixType Htilde = Fc0 * Harray[q].transpose();
        M(0, q)                 = Htilde.col(0).dot(g0);
        M(1, q)                 = Htilde.col(1).dot(g1);
        if constexpr (myDim == 2)
          M(2, q) = Htilde.col(1).dot(g0) + Htilde.col(0).dot(g1);
        else {
          const auto g2 = F.row(2);
          M(2, q)       = Htilde.col(2).dot(g2);
          M(3, q)       = Htilde.col(2).dot(g1) + Htilde.col(1).dot(g2);
          M(4, q)       = Htilde.col(2).dot(g0) + Htilde.col(0).dot(g2);
          M(5, q)       = Htilde.col(1).dot(g0) + Htilde.col(0).dot(g1);
        }
      }
      return M.eval();
    } else
      static_assert(Dune::AlwaysFalse<GEO>::value,
                    "firstDerivative can only be called with wrtCoeff as 0 and 1 indicating first derivative of the "
                    "Green-Lagrange strain w.r.t d and "
                    "alpha, respectively.");
  }

  /**
   * \brief Compute the second derivative of the Green-Lagrange strain w.r.t d or alpha for a given node and integration
   * point.
   *
   * \param geo The geometry object of the finite element.
   * \param uFunction The function representing the displacement field.
   * \param localBasis The local basis object.
   * \param gpPos The position of the integration point.
   * \param gpIndex The index of the integration point.
   * \param I The FE node index I (defaults to sNaN).
   * \param J The FE node index J (defaults to sNaN).
   * \param S The PK2 stress (in Voigt notation).
   * \param easFunction The EAS function.
   * \param alpha The coefficients of the EAS function.
   *
   * \tparam wrtCoeff An integer indicating the coefficient w.r.t which the second derivative of Green-Lagrange strain
   * is to be returned (0 = E,dd; 1 = E,aa; 2 = E,ad = E,da).
   * \tparam GEO The type of the geometry object.
   * \tparam ST The underlying scalar type.
   * \tparam EAST The type of the EAS function.
   *
   * \return The second derivative of the Green-Lagrange strain w.r.t d, alpha or the mixed derivative
   * (0 = E,dd; 1 = E,aa; 2 = E,ad = E,da).
   */
  template <int wrtCoeff, typename GEO, typename ST, typename EAST>
  static auto secondDerivative(const GEO& geo, const auto& uFunction, const auto& localBasis, const auto& gpIndex,
                               const Dune::FieldVector<double, GEO::mydimension>& gpPos,
                               const Eigen::Vector<ST, GEO::mydimension*(GEO::mydimension + 1) / 2>& S,
                               const EAST& easFunction, const auto& alpha, const int I = sNaN, const int J = sNaN) {
    constexpr int myDim              = GEO::mydimension;
    constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    using MatrixType     = Eigen::Matrix<ST, myDim, myDim>;
    const auto centerPos = centerPosition(geo);
    static_assert(myDim == 2 or myDim == 3,
                  "An enhancement of displacement gradient (transposed) can only be performed for the 2D or 3D case.");
    if constexpr (wrtCoeff == 0) { // E,dd
      const typename EAST::HType Harray = easFunction(gpPos);
      const auto gradUdI  = uFunction.evaluateDerivative(gpPos, wrt(spatialAll, coeff(I)), on(gridElement));
      const auto gradUdI0 = uFunction.evaluateDerivative(centerPos, wrt(spatialAll, coeff(I)), on(gridElement));
      const auto gradUdJ  = uFunction.evaluateDerivative(gpPos, wrt(spatialAll, coeff(J)), on(gridElement));
      const auto gradUdJ0 = uFunction.evaluateDerivative(centerPos, wrt(spatialAll, coeff(J)), on(gridElement));

      Eigen::Vector<ST, myDim> dNI, dNI0;
      Eigen::Vector<ST, myDim> dNJ, dNJ0;
      for (const auto i : Dune::range(myDim)) {
        dNI[i]  = getDiagonalEntry(gradUdI[i], 0);
        dNI0[i] = getDiagonalEntry(gradUdI0[i], 0);
        dNJ[i]  = getDiagonalEntry(gradUdJ[i], 0);
        dNJ0[i] = getDiagonalEntry(gradUdJ0[i], 0);
      }

      Eigen::Vector<ST, myDim> dNItilde = Eigen::Vector<ST, myDim>::Zero();
      Eigen::Vector<ST, myDim> dNJtilde = Eigen::Vector<ST, myDim>::Zero();

      for (const auto t : Dune::range(enhancedStrainSize)) {
        dNItilde += Harray[t] * (alpha[t] * dNI0);
        dNJtilde += Harray[t] * (alpha[t] * dNJ0);
      }

      dNI += dNItilde;
      dNJ += dNJtilde;
      ST val{0.0};
      if constexpr (myDim == 2)
        val = S[0] * dNI[0] * dNJ[0] + S[1] * dNI[1] * dNJ[1] + S[2] * (dNI[0] * dNJ[1] + dNI[1] * dNJ[0]);
      else
        val = S[0] * dNI[0] * dNJ[0] + S[1] * dNI[1] * dNJ[1] + S[2] * dNI[2] * dNJ[2] +
              S[3] * (dNI[2] * dNJ[1] + dNI[1] * dNJ[2]) + S[4] * (dNI[0] * dNJ[2] + dNI[2] * dNJ[0]) +
              S[5] * (dNI[0] * dNJ[1] + dNI[1] * dNJ[0]);

      const MatrixType kg = createScaledIdentityMatrix<ST, myDim, myDim>(val);
      return kg;
    } else if constexpr (wrtCoeff == 1) { // E,aa
      constexpr int enhancedStrainSize  = EAST::enhancedStrainSize;
      const typename EAST::HType Harray = easFunction(gpPos);
      auto kg                           = Eigen::Matrix<ST, enhancedStrainSize, enhancedStrainSize>::Zero().eval();
      const MatrixType Fc0 = MatrixType::Identity() + compatibleDisplacementGradient<GEO>(uFunction, centerPos);

      for (const auto i : Dune::range(myDim)) {
        for (const auto j : Dune::range(myDim)) {
          const auto sij = 0.5 * S[toVoigt<myDim>(i, j)];
          for (const auto p : Dune::range(enhancedStrainSize)) {
            const auto Hp = Fc0 * Harray[p].transpose();
            for (const auto q : Dune::range(enhancedStrainSize)) {
              const auto Hq = Fc0 * Harray[q].transpose();
              kg(p, q) += sij * (Hp.col(i).transpose() * Hq.col(j) + Hp.col(j).transpose() * Hq.col(i)).sum();
            }
          }
        }
      }
      return kg;
    } else if constexpr (wrtCoeff == 2) { // E,ad
      constexpr int enhancedStrainSize  = EAST::enhancedStrainSize;
      const typename EAST::HType Harray = easFunction(gpPos);
      auto kg                           = Eigen::Matrix<ST, enhancedStrainSize, myDim>::Zero().eval();
      const auto gradUdI   = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(I)), on(gridElement));
      const auto gradUdI0  = uFunction.evaluateDerivative(centerPos, wrt(spatialAll, coeff(I)), on(gridElement));
      const MatrixType Fc0 = MatrixType::Identity() + compatibleDisplacementGradient<GEO>(uFunction, centerPos);
      const MatrixType F   = computeDeformationGradient(geo, uFunction, gpPos, easFunction, alpha);

      Eigen::Vector<ST, myDim> dN0;
      for (const auto i : Dune::range(myDim))
        dN0[i] = getDiagonalEntry(gradUdI0[i], 0);

      Eigen::Vector<ST, myDim> dNtilde = Eigen::Vector<ST, myDim>::Zero();
      for (const auto t : Dune::range(enhancedStrainSize))
        dNtilde += Harray[t] * (alpha[t] * dN0);

      std::array<MatrixType, myDim> dNIdT, dNIdT0;
      for (const auto i : Dune::range(myDim)) {
        dNIdT[i].setZero();
        dNIdT0[i].setZero();
        for (const auto j : Dune::range(myDim)) {
          dNIdT[i](i, j)  = getDiagonalEntry(gradUdI[j], 0) + dNtilde[j];
          dNIdT0[i](i, j) = getDiagonalEntry(gradUdI0[j], 0);
        }
      }

      std::array<std::array<MatrixType, enhancedStrainSize>, myDim> dNXHtilde;
      for (auto& row : dNXHtilde)
        std::ranges::fill(row, MatrixType::Zero());
      for (const auto p : Dune::range(enhancedStrainSize)) {
        const auto Htilde = Harray[p].transpose();
        for (const auto q : Dune::range(myDim))
          dNXHtilde[q][p] += dNIdT0[q] * Htilde;
      }

      for (const auto i : Dune::range(myDim)) {
        for (const auto j : Dune::range(myDim)) {
          const auto sij = 0.5 * S[toVoigt<myDim>(i, j)];
          for (const auto p : Dune::range(enhancedStrainSize)) {
            const auto Hp = Fc0 * Harray[p].transpose();
            for (const auto q : Dune::range(myDim)) {
              kg(p, q) +=
                  sij * (Hp.col(i).transpose() * dNIdT[q].col(j) + Hp.col(j).transpose() * dNIdT[q].col(i) +
                         dNXHtilde[q][p].col(i).transpose() * F.col(j) + dNXHtilde[q][p].col(j).transpose() * F.col(i))
                            .sum();
            }
          }
        }
      }
      return kg;
    } else
      static_assert(
          Dune::AlwaysFalse<GEO>::value,
          "secondDerivative can only be called with wrtCoeff as 0, 1 and 2 indicating second derivative of the "
          "Green-Lagrange strain w.r.t d and "
          "alpha and the mixed derivative, respectively.");
  }

  /** \brief The name of the strain measure enhanced w.r.t EAS method. */
  static constexpr auto name() { return std::string("Displacement Gradient (Transposed)"); }

private:
  static constexpr int sNaN = std::numeric_limits<int>::signaling_NaN();

  template <typename GEO>
  static Dune::FieldVector<double, GEO::mydimension> centerPosition(const GEO& geo) {
    const auto& referenceElement = Dune::ReferenceElements<double, GEO::mydimension>::general(geo.type());
    const auto centerPos         = referenceElement.position(0, 0);
    return centerPos;
  }

  template <typename GEO>
  static auto compatibleDisplacementGradient(const auto& uFunction,
                                             const Dune::FieldVector<double, GEO::mydimension>& gpPos) {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    return toEigen(uFunction.evaluateDerivative(gpPos, wrt(spatialAll), on(gridElement))).eval();
  }

  template <typename GEO, typename EAST>
  static auto enhancedDisplacementGradient(const GEO& geo, const auto& uFunction,
                                           const Dune::FieldVector<double, GEO::mydimension>& gpPos,
                                           const EAST& easFunction, const auto& alpha) {
    using ST                          = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    constexpr int myDim               = GEO::mydimension;
    using MatrixType                  = Eigen::Matrix<ST, myDim, myDim>;
    const auto centerPos              = centerPosition(geo);
    const MatrixType Hc0              = compatibleDisplacementGradient<GEO>(uFunction, centerPos);
    const MatrixType Fc0              = (MatrixType::Identity() + Hc0).eval();
    const typename EAST::HType Harray = easFunction(gpPos);
    constexpr int enhancedStrainSize  = EAST::enhancedStrainSize;
    typename EAST::AnsatzType Htilde;
    Htilde.setZero(); // zeros returned if Harray is empty
    for (const auto p : Dune::range(enhancedStrainSize))
      Htilde += Harray[p] * alpha[p];
    const MatrixType Henhanced = (Fc0 * Htilde.transpose()).eval();
    return Henhanced;
  }

  template <typename GEO, typename EAST>
  static auto computeDisplacementGradient(const GEO& geo, const auto& uFunction,
                                          const Dune::FieldVector<double, GEO::mydimension>& gpPos,
                                          const EAST& easFunction, const auto& alpha) {
    using ST                = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    constexpr int myDim     = GEO::mydimension;
    using MatrixType        = Eigen::Matrix<ST, myDim, myDim>;
    const MatrixType Hc     = compatibleDisplacementGradient<GEO>(uFunction, gpPos);
    const MatrixType Htilde = enhancedDisplacementGradient(geo, uFunction, gpPos, easFunction, alpha);
    const MatrixType H      = Hc + Htilde;
    return H;
  }

  template <typename GEO, typename EAST>
  static auto computeDeformationGradient(const GEO& geo, const auto& uFunction,
                                         const Dune::FieldVector<double, GEO::mydimension>& gpPos,
                                         const EAST& easFunction, const auto& alpha) {
    using ST            = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    constexpr int myDim = GEO::mydimension;
    using MatrixType    = Eigen::Matrix<ST, myDim, myDim>;
    const MatrixType H  = computeDisplacementGradient(geo, uFunction, gpPos, easFunction, alpha);
    const MatrixType F  = (MatrixType::Identity() + H).eval();
    return F;
  }
};

} // namespace Ikarus::EAS

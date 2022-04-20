//
// Created by Alex on 20.04.2022.
//

#pragma once
#include "localFunctionInterface.hh"


namespace Ikarus {

template <typename Expr>
class LocalFunctionExpression : public LocalFunctionInterface<LocalFunctionExpression<Expr>> {
  using Base = LocalFunctionInterface<LocalFunctionExpression<Expr>>;

 public:
  friend Base;
  friend Expr;

  using Traits = LocalFunctionTraits<LocalFunctionExpression>;

  /** \brief Type used for coordinates */
  using ctype = typename Traits::ctype;
  /** \brief Dimension of the coeffs */
  static constexpr int valueSize = Traits::valueSize;
  /** \brief Dimension of the grid */
  static constexpr int gridDim = Traits::gridDim;
  /** \brief Type for the return value */
  using FunctionReturnType = typename Traits::FunctionReturnType;
  /** \brief Type for the directional derivatives */
  using AlongType = typename Traits::AlongType;
  /** \brief Type for the coordinates to store the return value */
  using GlobalE = typename FunctionReturnType::CoordinateType;
  /** \brief Type for the Jacobian matrix */
  using Jacobian = typename Traits::Jacobian;
  /** \brief Type for a column of the Jacobian matrix */
  using JacobianColType = typename Traits::JacobianColType;
  /** \brief Type for the derivatives wrt. the coeffiecients */
  using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
  /** \brief Type for ansatz function values */
  using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
  /** \brief Type for the Jacobian of the ansatz function values */
  using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

  private:
  const Expr& underlying() const {return static_cast<Expr const&>(*this);}


  FunctionReturnType evaluateFunctionImpl(const Eigen::VectorXd& N) const { return underlying().evaluateFunctionImpl(N); }

  Jacobian evaluateDerivativeWRTSpaceAllImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN) const { return underlying().evaluateDerivativeWRTSpaceAllImpl(N,dN); }

  JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN,
                                                       int spaceIndex) const { return underlying().evaluateDerivativeWRTSpaceSingleImpl(N,dN,spaceIndex); }

  CoeffDerivMatrix evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N,
                                                   const AnsatzFunctionJacobian& dN,
                                                   int coeffsIndex) const { return underlying().evaluateDerivativeWRTCoeffsImpl(N,dN,coeffsIndex); }

  CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffs(const AnsatzFunctionType& N,
                                                     const AnsatzFunctionJacobian& dN,
                                                     const AlongType& along,
                                                     const std::array<size_t, 2>& coeffsIndex) const { return underlying().evaluateSecondDerivativeWRTCoeffs(N,dN,along,coeffsIndex); }

  std::array<CoeffDerivMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
      const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN, int coeffsIndex) const { return underlying().evaluateDerivativeWRTCoeffsANDSpatialImpl(N,dN,coeffsIndex); }

  CoeffDerivMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                                   const AnsatzFunctionJacobian& dN,
                                                                   int coeffsIndex, const int spatialIndex) const { return underlying().evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(N,dN,coeffsIndex,spatialIndex); }

  std::array<CoeffDerivMatrix, gridDim> evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
      const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN, const AlongType& along,
      const std::array<size_t, 2>& coeffsIndex) const  { return underlying().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(N,dN,along,coeffsIndex); }

  CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
      const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN, const AlongType& along,
      const std::array<size_t, 2>& coeffsIndex, const int spatialIndex) const { return underlying().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(N,dN,along,coeffsIndex,spatialIndex); }

  const auto& basis() const
  {
    return underlying().basis();
  }
};

template <typename Expr>
struct LocalFunctionTraits<LocalFunctionExpression<Expr>> : public LocalFunctionTraits<Expr>{
  using Base =  LocalFunctionTraits<Expr>;

};

}  // namespace Ikarus

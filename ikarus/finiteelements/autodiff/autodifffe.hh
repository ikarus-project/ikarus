// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file autoDiffFE.hh
 * \brief Contains the AutoDiffFE class, an automatic differentiation wrapper for finite elements.
 */

#pragma once

#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/fetraits.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

/**
 * \brief AutoDiffFE class, an automatic differentiation wrapper for finite elements.
 *
 * \tparam FEImpl The type of the original finite element, which does not implement the derivatives
 * \tparam FER Type of the Finite Element Requirements.
 * \tparam useEigenRef A boolean indicating whether to use Eigen::Ref for references types in calculateMatrix,...
 * \tparam forceAutoDiff A boolean indicating whether to force the use of automatic differentiation, even when the
 * real element implements the derivatives.
 */
template <typename FEImpl, typename FER = FERequirements<>, bool useEigenRef = false, bool forceAutoDiff = false>
class AutoDiffFE : public FEImpl
{
public:
  using RealFE            = FEImpl;                             ///< Type of the base finite element.
  using Basis             = typename RealFE::Basis;             ///< Type of the basis.
  using Traits            = FETraits<Basis, FER, useEigenRef>;  ///< Type traits for local view.
  using LocalView         = typename Traits::LocalView;         ///< Type of the local view.
  using Element           = typename Traits::Element;           ///< Type of the element.
  using FERequirementType = typename Traits::FERequirementType; ///< Type of the Finite Element Requirements.

  /**
   * \brief Calculate the matrix associated with the finite element.
   *
   * \param req Finite Element Requirements.
   * \param h Matrix to be calculated.
   */
  void calculateMatrix(const FERequirementType& req, typename Traits::template MatrixType<> h) const {
    // real element implements calculateMatrix by itself, then we simply forward the call
    if constexpr (requires { RealFE::calculateMatrix(req, h); } and not forceAutoDiff) {
      RealFE::calculateMatrix(req, h);
    } else if constexpr (requires {
                           this->template calculateVectorImpl<autodiff::dual>(
                               req, std::declval<typename Traits::template VectorType<autodiff::dual>>(),
                               std::declval<const Eigen::VectorXdual&>());
                         }) {
      // real element implements calculateVector by itself, therefore we only need first order derivatives
      Eigen::VectorXdual dx(this->localView().size());
      Eigen::VectorXdual g(this->localView().size());
      dx.setZero();
      auto f = [&](auto& x) -> auto& {
        g.setZero();
        this->template calculateVectorImpl<autodiff::dual>(req, g, x);
        return g;
      };
      jacobian(f, autodiff::wrt(dx), at(dx), g, h);
    } else if constexpr (requires {
                           this->template calculateScalarImpl<autodiff::dual2nd>(
                               req, std::declval<typename Traits::template VectorType<autodiff::dual2nd>>());
                         }) {
      // real element implements calculateScalar by itself, therefore we need second order derivatives
      Eigen::VectorXdual2nd dx(this->localView().size());
      Eigen::VectorXd g;
      autodiff::dual2nd e;
      dx.setZero();
      auto f = [&](auto& x) { return this->template calculateScalarImpl<autodiff::dual2nd>(req, x); };
      hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
    } else
      static_assert(Dune::AlwaysFalse<AutoDiffFE>::value,
                    "Appropriate calculateScalarImpl or calculateVectorImpl functions are not implemented for the "
                    "chosen element.");
  }

  /**
   * \brief Calculate the vector associated with the finite element.
   *
   * \param req Finite Element Requirements.
   * \param g Vector to be calculated.
   */
  void calculateVector(const FERequirementType& req, typename Traits::template VectorType<> g) const {
    // real element implements calculateVector by itself, then we simply forward the call
    if constexpr (requires {
                    this->template calculateVectorImpl<double>(
                        req, std::declval<typename Traits::template VectorType<double>>(),
                        std::declval<const Eigen::VectorXd&>());
                  } and not forceAutoDiff) {
      return this->template calculateVectorImpl<double>(req, g);
    } else if constexpr (requires {
                           this->template calculateScalarImpl<autodiff::dual>(
                               req, std::declval<const Eigen::VectorXdual&>());
                         }) {
      // real element implements calculateScalar by itself, therefore we need first order derivatives
      Eigen::VectorXdual dx(this->localView().size());
      dx.setZero();
      autodiff::dual e;
      auto f = [&](auto& x) { return this->template calculateScalarImpl<autodiff::dual>(req, x); };
      gradient(f, autodiff::wrt(dx), at(dx), e, g);
    } else
      static_assert(Dune::AlwaysFalse<AutoDiffFE>::value,
                    "Appropriate calculateScalarImpl function is not implemented for the "
                    "chosen element.");
  }

  /**
   * \brief Calculate the local system associated with the finite element.
   *
   * \param req Finite Element Requirements.
   * \param h Matrix to be calculated.
   * \param g Vector to be calculated.
   */
  void calculateLocalSystem(const FERequirementType& req, typename Traits::template MatrixType<> h,
                            typename Traits::template VectorType<> g) const {
    Eigen::VectorXdual2nd dx(this->localView().size());
    dx.setZero();
    auto f = [&](auto& x) { return this->calculateScalarImpl(req, x); };
    hessian(f, autodiff::wrt(dx), at(dx), g, h);
  }

  /**
   * \brief Calculate the scalar value associated with the finite element.
   *
   * \param par Finite Element Requirements.
   * \return The calculated scalar value.
   */
  [[nodiscard]] double calculateScalar(const FERequirementType& par) const {
    // real element implements calculateScalar by itself, then we simply forward the call
    if constexpr (requires { RealFE::calculateScalar(par); }) {
      return RealFE::calculateScalar(par);
    } else if constexpr (requires { this->calculateScalarImpl(par); }) {
      // real element only implements the protected calculateScalarImpl by itself, thus we call that one.
      return this->calculateScalarImpl(par);
    } else {
      static_assert(Dune::AlwaysFalse<AutoDiffFE>::value,
                    "Appropriate calculateScalar and calculateScalarImpl functions are not implemented for the "
                    "chosen element.");
    }
  }

  /**
   * \brief Get the reference to the base finite element.
   *
   * \return The reference to the base finite element.
   */
  const RealFE& realFE() const { return *this; }

  /**
   * \brief Constructor for the AutoDiffFE class.
   * Forward the construction to the underlying element
   *
   * \tparam Args Variadic template for constructor arguments.
   * \param args Constructor arguments.
   */
  template <typename... Args>
  explicit AutoDiffFE(Args&&... args)
      : RealFE{std::forward<Args>(args)...} {}
};
} // namespace Ikarus

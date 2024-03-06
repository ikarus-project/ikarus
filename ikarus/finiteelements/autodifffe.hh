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
#include <ikarus/utils/traits.hh>

namespace Ikarus {

/**
 * \brief AutoDiffFE class, an automatic differentiation wrapper for finite elements.
 *
 * \tparam FEImpl The type of the original finite element, which does not implement the derivatives
 * \tparam forceAutoDiff A boolean indicating whether to force the use of automatic differentiation, even when the
 * real element implements the derivatives.
 */
template <typename FEImpl, bool forceAutoDiff = false>
class AutoDiffFE : public FEImpl
{
public:
  using RealFE            = FEImpl;                             ///< Type of the base finite element.
  using BasisHandler      = typename RealFE::BasisHandler;      ///< Type of the basis handler.
  using Traits            = typename RealFE::Traits;            ///< Type traits for local view.
  using LocalView         = typename Traits::LocalView;         ///< Type of the local view.
  using Element           = typename Traits::Element;           ///< Type of the element.
  using FERequirementType = typename Traits::FERequirementType; ///< Type of the Finite Element Requirements.
private:
  using Mixin = FEImpl::Mixin;

public:
  /**
   * \brief Calculate the matrix associated with the finite element.
   *
   * \param self The finite element based on AutoDiff itself.
   * \param par Finite Element Requirements.
   * \param h Matrix to be calculated.
   */
  friend void calculateMatrix(const AutoDiffFE& self, const FERequirementType& par,
                              typename Traits::template MatrixType<> h) {
    self.calculateMatrix(par, h);
  }

  /**
   * \brief Calculate the vector associated with the finite element.
   *
   * \param self The finite element based on AutoDiff itself.
   * \param par Finite Element Requirements.
   * \param g Vector to be calculated.
   */
  friend void calculateVector(const AutoDiffFE& self, const FERequirementType& par,
                              typename Traits::template VectorType<double> g) {
    self.calculateVector(par, g);
  }

  /**
   * \brief Calculate the local system associated with the finite element.
   *
   * \param self The finite element based on AutoDiff itself.
   * \param par Finite Element Requirements.
   * \param h Matrix to be calculated.
   * \param g Vector to be calculated.
   */
  friend void calculateLocalSystem(const AutoDiffFE& self, const FERequirementType& par,
                                   typename Traits::template MatrixType<> h, typename Traits::template VectorType<> g) {
    self.calculateLocalSystem(par, h, g);
  }

  /**
   * \brief Calculate the scalar value associated with the finite element.
   *
   * \param self The finite element based on AutoDiff itself.
   * \param par Finite Element Requirements.
   * \return The calculated scalar value.
   */
  friend auto calculateScalar(const AutoDiffFE& self, const FERequirementType& par) {
    return self.calculateScalar(par);
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

private:
  void calculateMatrix(const FERequirementType& req, typename Traits::template MatrixType<> h) const {
    // real element implements calculateMatrix by itself, then we simply forward the call

    if constexpr (requires(Eigen::VectorXd v) {
                    static_cast<const Mixin&>(std::declval<AutoDiffFE>())
                        .template calculateMatrixImpl<double>(req, h, v);
                  } and not forceAutoDiff) {
      return Mixin::template calculateMatrixImpl<double>(req, h);
    } else if constexpr (requires(Eigen::VectorXdual v) {
                           static_cast<const Mixin&>(std::declval<AutoDiffFE>())
                               .template calculateVectorImpl<autodiff::dual>(
                                   req, std::declval<typename Traits::template VectorType<autodiff::dual>>(), v);
                         }) {
      // real element implements calculateVector by itself, therefore we only need first order derivatives
      Eigen::VectorXdual dx(this->localView().size());
      Eigen::VectorXdual g(this->localView().size());
      dx.setZero();
      auto f = [this, &req, &g](auto& x) -> auto& {
        // Since req is const as a function argument, we can not make this lambda capture by mutable reference
        // But we have to do this since for efficiency reason we reuse the g vector
        // therefore, the only remaining option is to cast the const away from g
        Eigen::VectorXdual& gref = const_cast<Eigen::VectorXdual&>(g);
        gref.setZero();
        Mixin::template calculateVectorImpl<autodiff::dual>(req, gref, x);
        return g;
      };
      jacobian(f, autodiff::wrt(dx), at(dx), g, h);
    } else if constexpr (requires(typename Traits::template VectorType<autodiff::dual2nd> v) {
                           static_cast<const Mixin&>(std::declval<AutoDiffFE>())
                               .template calculateScalarImpl<autodiff::dual2nd>(req, v);
                         }) {
      // real element implements calculateScalar by itself, therefore we need second order derivatives
      Eigen::VectorXdual2nd dx(this->localView().size());
      Eigen::VectorXd g;
      autodiff::dual2nd e;
      dx.setZero();
      auto f = [this, &req](auto& x) { return Mixin::template calculateScalarImpl<autodiff::dual2nd>(req, x); };
      hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
    } else
      static_assert(Dune::AlwaysFalse<AutoDiffFE>::value,
                    "Appropriate calculateScalarImpl or calculateVectorImpl functions are not implemented for the "
                    "chosen element.");
  }

  void calculateVector(const FERequirementType& req, typename Traits::template VectorType<> g) const {
    // real element implements calculateVector by itself, then we simply forward the call
    if constexpr (requires {
                    static_cast<const Mixin&>(std::declval<AutoDiffFE>())
                        .template calculateVectorImpl<double>(
                            req, std::declval<typename Traits::template VectorType<double>>(),
                            std::declval<const Eigen::VectorXd&>());
                  } and not forceAutoDiff) {
      return Mixin::template calculateVectorImpl<double>(req, g);
    } else if constexpr (requires {
                           static_cast<const Mixin&>(std::declval<AutoDiffFE>())
                               .template calculateScalarImpl<autodiff::dual>(req,
                                                                             std::declval<const Eigen::VectorXdual&>());
                         }) {
      // real element implements calculateScalar by itself but no calculateVectorImpl, therefore we need first order
      // derivatives
      Eigen::VectorXdual dx(this->localView().size());
      dx.setZero();
      autodiff::dual e;
      auto f = [this, &req](auto& x) { return Mixin::template calculateScalarImpl<autodiff::dual>(req, x); };
      gradient(f, autodiff::wrt(dx), at(dx), e, g);
    } else
      static_assert(Dune::AlwaysFalse<AutoDiffFE>::value,
                    "Appropriate calculateScalarImpl function is not implemented for the "
                    "chosen element.");
  }

  [[nodiscard]] double calculateScalar(const FERequirementType& par) const {
    // real element implements calculateScalar by itself, then we simply forward the call
    if constexpr (requires {
                    static_cast<const Mixin&>(std::declval<AutoDiffFE>()).template calculateScalarImpl<double>(par);
                  }) {
      Mixin::template calculateScalarImpl<double>(par);
      // real element only implements the protected calculateScalarImpl by itself, thus we call that one.
      return Mixin::template calculateScalarImpl<double>(par);
    } else {
      static_assert(Dune::AlwaysFalse<AutoDiffFE>::value,
                    "Appropriate calculateScalar and calculateScalarImpl functions are not implemented for the "
                    "chosen element.");
    }
  }

  void calculateLocalSystem(const FERequirementType& req, typename Traits::template MatrixType<> h,
                            typename Traits::template VectorType<> g) const {
    Eigen::VectorXdual2nd dx(this->localView().size());
    dx.setZero();
    auto f = [&](auto& x) { return Mixin::calculateScalarImpl(req, x); };
    hessian(f, autodiff::wrt(dx), at(dx), g, h);
  }
};
} // namespace Ikarus

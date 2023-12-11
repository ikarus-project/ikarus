// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

////
//
#pragma once
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {
  template <typename RealElement, typename FERequirementType_ = FErequirements<>, bool useEigenRef = false,
            bool forceAutoDiff = false>
  class AutoDiffFE : public RealElement {
  public:
    using Base              = RealElement;
    using Basis             = Base::Basis;
    using LocalView         = typename Basis::FlatBasis::LocalView;
    using Traits            = TraitsFromLocalView<LocalView, useEigenRef>;
    using Element           = typename LocalView::Element;
    using FERequirementType = FERequirementType_;

    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> h) const {
      if constexpr (requires { RealElement::calculateMatrix(par, h); } and not forceAutoDiff) {
        RealElement::calculateMatrix(par, h);
      } else if constexpr (requires {
                             this->template calculateVectorImpl<autodiff::dual>(
                                 par, std::declval<typename Traits::template VectorType<autodiff::dual>>(),
                                 std::declval<const Eigen::VectorXdual&>());
                           }) {
        /// This is only valid if the external forces are independent of displacements, for e.g., no follower forces are
        /// applied
        Eigen::VectorXdual dx(this->localView().size());
        Eigen::VectorXdual g(this->localView().size());
        dx.setZero();
        auto f = [&](auto& x) -> auto& {
          g.setZero();
          this->template calculateVectorImpl<autodiff::dual>(par, g, x);
          return g;
        };
        jacobian(f, autodiff::wrt(dx), at(dx), g, h);
      } else if constexpr (requires {
                             this->template calculateScalarImpl<autodiff::dual2nd>(
                                 par, std::declval<typename Traits::template VectorType<autodiff::dual2nd>>());
                           }) {
        Eigen::VectorXdual2nd dx(this->localView().size());
        Eigen::VectorXd g;
        autodiff::dual2nd e;
        dx.setZero();
        auto f = [&](auto& x) { return this->template calculateScalarImpl<autodiff::dual2nd>(par, x); };
        hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
      } else
        static_assert(Ikarus::Std::DummyFalse<AutoDiffFE>::value,
                      "Appropriate calculateScalarImpl or calculateVectorImpl functions are not implemented for the "
                      "chosen element.");
    }

    inline void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> g) const {
      if constexpr (requires {
                      this->template calculateVectorImpl<double>(
                          par, std::declval<typename Traits::template VectorType<double>>(),
                          std::declval<const Eigen::VectorXd&>());
                    }
                    and not forceAutoDiff) {
        return this->template calculateVectorImpl<double>(par, g);
      } else if constexpr (requires {
                             this->template calculateScalarImpl<autodiff::dual>(
                                 par, std::declval<const Eigen::VectorXdual&>());
                           }) {
        Eigen::VectorXdual dx(this->localView().size());
        dx.setZero();
        autodiff::dual e;
        auto f = [&](auto& x) { return this->template calculateScalarImpl<autodiff::dual>(par, x); };
        gradient(f, autodiff::wrt(dx), at(dx), e, g);
      } else
        static_assert(Ikarus::Std::DummyFalse<AutoDiffFE>::value,
                      "Appropriate calculateScalarImpl function is not implemented for the "
                      "chosen element.");
    }

    void calculateLocalSystem(const FERequirementType& par, typename Traits::template MatrixType<> h,
                              typename Traits::template VectorType<> g) const {
      Eigen::VectorXdual2nd dx(this->localView().size());
      dx.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), g, h);
    }

    [[nodiscard]] double calculateScalar(const FERequirementType& par) const {
      if constexpr (requires { RealElement::calculateScalar(par); }) {
        return RealElement::calculateScalar(par);
      } else if constexpr (requires { this->calculateScalarImpl(par); }) {
        return this->calculateScalarImpl(par);
      } else
        static_assert(Ikarus::Std::DummyFalse<AutoDiffFE>::value,
                      "Appropriate calculateScalar and calculateScalarImpl functions are not implemented for the "
                      "chosen element.");
    }

    const RealElement& getFE() const { return *this; }

    template <typename... Args>
    explicit AutoDiffFE(Args&&... args) : RealElement{std::forward<Args>(args)...} {}
  };
}  // namespace Ikarus

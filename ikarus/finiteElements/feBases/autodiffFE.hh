// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

////
//
#pragma once
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/physicsHelper.hh>

namespace Ikarus {
  template <typename RealElement, typename FERequirementType_ = FErequirements<>, bool useEigenRef = false>
  class AutoDiffFE : public RealElement {
  public:
    using Base              = RealElement;
    using Basis             = Base::Basis;
    using LocalView         = typename Basis::FlatBasis::LocalView;
    using Traits            = TraitsFromLocalView<LocalView, useEigenRef>;
    using Element           = typename LocalView::Element;
    using FERequirementType = FERequirementType_;

    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<>& h) const {
      if constexpr (requires {
                      this->calculateVectorImpl(par, std::declval<const Eigen::VectorXdual&>(),
                                                std::declval<typename Traits::template VectorType<autodiff::dual>>());
                    }) {
        /// This is only valid if the external forces are independent of displacements, for e.g., no follower forces are
        /// applied
        Eigen::VectorXdual dx(this->localView().size());
        Eigen::VectorXdual g(this->localView().size());
        dx.setZero();
        auto f = [&](auto& x) -> auto& {
          g.setZero();
          this->calculateVectorImpl(par, x, g);
          return g;
        };
        jacobian(f, autodiff::wrt(dx), at(dx), g, h);
      } else if constexpr (requires {
                             this->calculateScalarImpl(
                                 par, std::declval<typename Traits::template VectorType<autodiff::dual2nd>>());
                           }) {
        Eigen::VectorXdual2nd dx(this->localView().size());
        Eigen::VectorXd g;
        autodiff::dual2nd e;
        dx.setZero();
        auto f = [&](auto& x) { return this->template calculateScalarImpl<autodiff::dual2nd>(par, x); };
        hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
      } else
        DUNE_THROW(Dune::NotImplemented,
                   "Appropriate calculateScalarImpl and calculateVectorImpl functions are not implemented for the "
                   "chosen element.");
    }

    void calculateVector(const FERequirementType& par, typename Traits::template VectorType<>& g) const {
      if constexpr (requires {
                      this->calculateVectorImpl(par, std::declval<const Eigen::VectorXd&>(),
                                                std::declval<typename Traits::template VectorType<double>>());
                    }) {
        Eigen::VectorXd dx(this->localView().size());
        dx.setZero();
        return this->calculateVectorImpl(par, dx, g);
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
        DUNE_THROW(Dune::NotImplemented,
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
      } else if constexpr (requires { this->calculateScalarImpl(par, std::declval<const Eigen::VectorXd&>()); }) {
        Eigen::VectorXd dx(this->localView().size());
        dx.setZero();
        return this->calculateScalarImpl(par, dx);
      } else
        DUNE_THROW(Dune::NotImplemented,
                   "Appropriate calculateScalar and calculateScalarImpl functions are not implemented for the "
                   "chosen element.");
    }

    const RealElement& getFE() const { return *this; }

    template <typename... Args>
    explicit AutoDiffFE(Args&&... args) : RealElement{std::forward<Args>(args)...} {}

    //    explicit AutoDiffFE(RealElement& other) : RealElement(other) {
    //      if constexpr (requires { this->setEASType(int{}); }) {
    //        const auto& numberOfEASParameters = other.getNumberOfEASParameters();
    //        this->setEASType(numberOfEASParameters);
    //      }
    //    }
  };
}  // namespace Ikarus

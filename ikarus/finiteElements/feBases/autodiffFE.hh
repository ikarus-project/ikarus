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

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      if constexpr (requires {
                      this->calculateVectorImpl(par, typename Traits::VectorType{}, Eigen::VectorX<double>{});
                    }) {
        /// This is only valid if the external forces are independent of displacements, for e.g., no follower forces are
        /// applied
        std::cout << "Hello calculateMatrix from calculateVectorImpl\n";
        Eigen::VectorXdual dx(this->localView().size());
        typename Traits::VectorType g(this->localView().size());
        dx.setZero();
        autodiff::dual e;
        auto f = [&](auto& x) { return this->calculateVectorImpl(par, g, x); };
        gradient(f, autodiff::wrt(dx), at(dx), e, h);
      } else if constexpr (requires { this->calculateScalarImpl(par, Eigen::VectorX<double>{}); }) {
        std::cout << "Hello calculateMatrix from calculateScalarImpl\n";
        Eigen::VectorXdual2nd dx(this->localView().size());
        Eigen::VectorXd g;
        autodiff::dual2nd e;
        dx.setZero();
        auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
        hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
      } else
        DUNE_THROW(Dune::NotImplemented,
                   "Appropriate calculateScalarImpl and calculateVectorImpl functions are not implemented for the "
                   "chosen element.");
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      if constexpr (requires { this->calculateScalarImpl(par, Eigen::VectorX<double>{}); }) {
        Eigen::VectorXdual dx(this->localView().size());
        dx.setZero();
        autodiff::dual e;
        auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
        gradient(f, autodiff::wrt(dx), at(dx), e, g);
      } else
        DUNE_THROW(Dune::NotImplemented,
                   "Appropriate calculateScalarImpl function is not implemented for the "
                   "chosen element.");
    }

    void calculateLocalSystem(const FERequirementType& par, typename Traits::MatrixType& h,
                              typename Traits::VectorType& g) const {
      Eigen::VectorXdual2nd dx(this->localView().size());
      dx.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), g, h);
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      if constexpr (requires { this->calculateScalarImpl(par, Eigen::VectorX<double>{}); }) {
        Eigen::VectorXd dx(this->localView().size());
        dx.setZero();
        return this->calculateScalarImpl(par, dx);
      } else
        DUNE_THROW(Dune::NotImplemented,
                   "Appropriate calculateScalarImpl function is not implemented for the "
                   "chosen element.");
    }

    template <typename... Args>
    explicit AutoDiffFE(Args&&... args) : RealElement{std::forward<Args>(args)...} {}
  };
}  // namespace Ikarus

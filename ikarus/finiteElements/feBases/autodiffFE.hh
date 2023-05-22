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
  /// This element can not be used on its own but it should be inherited from
  /// The class constructor can only be called from the templated class.
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
      Eigen::VectorXdual2nd dx(this->localView().size());
      Eigen::VectorXd g;
      autodiff::dual2nd e;
      dx.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      Eigen::VectorXdual dx(this->localView().size());
      dx.setZero();
      autodiff::dual e;
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      gradient(f, autodiff::wrt(dx), at(dx), e, g);
    }

    void calculateLocalSystem(const FERequirementType& par, typename Traits::MatrixType& h,
                              typename Traits::VectorType& g) const {
      Eigen::VectorXdual2nd dx(this->localView().size());
      dx.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), g, h);
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(this->localView().size());
      dx.setZero();
      return this->calculateScalarImpl(par, dx);
    }

    template <typename... Args>
    explicit AutoDiffFE(Args&&... args) : RealElement{std::forward<Args>(args)...} {}
  };
}  // namespace Ikarus

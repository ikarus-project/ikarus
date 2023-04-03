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
  ////This element can not be used on its own but it should be inherited from
  //// The class constructor can only be called from the templated class.
  template <typename RealElement, typename Basis, typename FERequirementTypeImpl = FErequirements<Eigen::VectorXd>>
  class AutoDiffFE {
  public:
    using LocalView             = typename Basis::LocalView;
    using Traits                = TraitsFromLocalView<LocalView>;
    using GridElementEntityType = typename LocalView::Element;
    using FERequirementType     = FErequirements<Eigen::VectorXd>;
    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      Eigen::VectorXdual2nd dx(localDofSize);
      Eigen::VectorXd g;
      autodiff::dual2nd e;
      dx.setZero();
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      Eigen::VectorXdual dx(localDofSize);
      dx.setZero();
      autodiff::dual e;
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      gradient(f, autodiff::wrt(dx), at(dx), e, g);
    }

    void calculateLocalSystem(const FERequirementType& par, typename Traits::MatrixType& h,
                              typename Traits::VectorType& g) const {
      Eigen::VectorXdual2nd dx(localDofSize);
      dx.setZero();
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), g, h);
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(localDofSize);
      dx.setZero();

      return this->underlying().calculateScalarImpl(par, dx);
    }

    [[nodiscard]] size_t size() const { return localDofSize; }
    const GridElementEntityType& getEntity() { return localView_.element(); }
    const LocalView& localView() const { return localView_; }
    LocalView& localView() { return localView_; }

  protected:
    explicit AutoDiffFE(const Basis& basis, const typename LocalView::Element& element)
        : localView_{basis.localView()} {
      localView_.bind(element);
      localDofSize = localView_.size();
    }

  private:
    RealElement const& underlying() const  // CRTP
    {
      return static_cast<RealElement const&>(*this);
    }
    LocalView localView_;
    int localDofSize{};
  };

}  // namespace Ikarus

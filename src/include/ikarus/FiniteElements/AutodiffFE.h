////
//// Created by Alex on 05.07.2021.
////
//
#pragma once
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <ikarus/FiniteElements/physicsHelper.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>

namespace Ikarus {
  ////This element can not be used on its own but it should be inherited from
  //// The class constructor can only be called from the templated class.
  template <typename RealElement>
  class AutoDiffFE {

  public:
    template <typename LocalView, typename Derived>
    static double calculateScalar(const LocalView& localView, const Eigen::MatrixBase<Derived>& d, const double& lambda=0) {
      using namespace autodiff;
      Eigen::Vector<double,Derived::RowsAtCompileTime> dx(localView.size());
      return RealElement::calculateScalarImpl(localView, d, dx,lambda);
    }

    template <typename LocalView, typename Derived>
    static auto calculateMatrix(const LocalView& localView, const Eigen::MatrixBase<Derived>& d, const double& lambda=0) {
      using namespace autodiff;
      Eigen::Vector<dual2nd,Derived::RowsAtCompileTime> dx(localView.size());
      dx.setZero();
      auto f = [&](auto& x) { return RealElement::calculateScalarImpl(localView, d, x,lambda); };
      return hessian(f, wrt(dx), at(dx));
    }

    template <typename LocalView, typename Derived>
    static auto calculateVector(const LocalView& localView, const Eigen::MatrixBase<Derived>& d, const double& lambda=0) {
      using namespace autodiff;
      Eigen::Vector<dual,Derived::RowsAtCompileTime> dx(localView.size());
      dx.setZero();
      auto f = [&](auto& x) { return RealElement::calculateScalarImpl(localView, d, x,lambda); };
      return gradient(f, wrt(dx), at(dx));
    }
  };


  template <typename RealElement,typename LocalView>
  class AutoDiffFEClean {

  public:
//    using LocalView = typename RealElement::LocalView;
    using Traits = TraitsFromLocalView<LocalView>;
    AutoDiffFEClean(const LocalView& localView)
        :  localView_{localView} {}
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    [[nodiscard]] typename Traits::MatrixType calculateMatrix(const FERequirementType &par) const {
      Eigen::VectorXdual2nd dx(localView_.size());
      dx.setZero();
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      return hessian(f, wrt(dx), at(dx));
    }

    [[nodiscard]] typename Traits::VectorType calculateVector(const FERequirementType &par) const {
      Eigen::VectorXdual dx(localView_.size());
      dx.setZero();
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      return gradient(f, wrt(dx), at(dx));
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType &par) const {
      Eigen::VectorXd dx(localView_.size());

      return this->underlying().calculateScalarImpl(par, dx);
    }

  private:
    RealElement const& underlying() const //CRTP
    {
            return  static_cast<RealElement const& >(*this);
    }

    LocalView localView_;
  };

}  // namespace Ikarus
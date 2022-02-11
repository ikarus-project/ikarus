////
//// Created by Alex on 05.07.2021.
////
//
#pragma once
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace Ikarus {
  ////This element can not be used on its own but it should be inherited from
  //// The class constructor can only be called from the templated class.
  template <typename RealElement>
  class AutoDiffFE {

  public:
    template <typename LocalView, typename Derived>
    static double calculateScalar(const LocalView& localView, const Eigen::MatrixBase<Derived>& d, const double& lambda) {
      using namespace autodiff;
      Eigen::Vector<double,Derived::RowsAtCompileTime> dx(localView.size());
      return RealElement::template calculateScalarImpl(localView, d, dx,lambda);
    }

    template <typename LocalView, typename Derived>
    static auto calculateMatrix(const LocalView& localView, const Eigen::MatrixBase<Derived>& d, const double& lambda) {
      using namespace autodiff;
      Eigen::Vector<dual2nd,Derived::RowsAtCompileTime> dx(localView.size());
      dx.setZero();
      auto f = [&](auto& x) { return RealElement::template calculateScalarImpl(localView, d, x,lambda); };
      return hessian(f, wrt(dx), at(dx));
    }

    template <typename LocalView, typename Derived>
    static auto calculateVector(const LocalView& localView, const Eigen::MatrixBase<Derived>& d, const double& lambda) {
      using namespace autodiff;
      Eigen::Vector<dual,Derived::RowsAtCompileTime> dx(localView.size());
      dx.setZero();
      auto f = [&](auto& x) { return RealElement::template calculateScalarImpl(localView, d, x,lambda); };
      return gradient(f, wrt(dx), at(dx));
    }
  };

}  // namespace Ikarus
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
    template <typename LocalView>
    static double calculateScalar(const LocalView& localView, const Eigen::VectorXd& d) {
      return RealElement::template calculateScalarImpl<double>(localView, d);
    }

    template <typename LocalView>
    static auto calculateMatrix(const LocalView& localView, const Eigen::VectorXd& d) {
      Eigen::Vector4dual2nd dx(localView.size());
      dx.setZero();
      auto f = [&](auto& x) { return RealElement::template calculateScalarImpl(localView, d, x); };
      return hessian(f, wrt(dx), at(dx));
    }

    template <typename LocalView>
    static auto calculateVector(const LocalView& localView, const Eigen::VectorXd& d) {
      Eigen::Vector4dual dx(localView.size());
      dx.setZero();
      auto f = [&](auto& x) { return RealElement::template calculateScalarImpl(localView, d, x); };
      return gradient(f, wrt(dx), at(dx));
    }
  };

}  // namespace Ikarus
////
//// Created by Alex on 05.07.2021.
////
//
#pragma once
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/physicsHelper.h>

namespace Ikarus {
  ////This element can not be used on its own but it should be inherited from
  //// The class constructor can only be called from the templated class.
  template <typename RealElement>
  class AutoDiffFE {
  public:
    template <typename LocalView, typename Derived>
    static double calculateScalar(const LocalView& localView, const Eigen::MatrixBase<Derived>& d,
                                  const double& lambda = 0) {
      using namespace autodiff;
      Eigen::Vector<double, Derived::RowsAtCompileTime> dx(localView.size());
      return RealElement::calculateScalarImpl(localView, d, dx, lambda);
    }

    template <typename LocalView, typename Derived>
    static auto calculateMatrix(const LocalView& localView, const Eigen::MatrixBase<Derived>& d,
                                const double& lambda = 0) {
      using namespace autodiff;
      Eigen::Vector<dual2nd, Derived::RowsAtCompileTime> dx(localView.size());
      dx.setZero();
      auto f = [&](auto& x) { return RealElement::calculateScalarImpl(localView, d, x, lambda); };
      return hessian(f, wrt(dx), at(dx));
    }

    template <typename LocalView, typename Derived>
    static auto calculateVector(const LocalView& localView, const Eigen::MatrixBase<Derived>& d,
                                const double& lambda = 0) {
      using namespace autodiff;
      Eigen::Vector<dual, Derived::RowsAtCompileTime> dx(localView.size());
      dx.setZero();
      auto f = [&](auto& x) { return RealElement::calculateScalarImpl(localView, d, x, lambda); };
      return gradient(f, wrt(dx), at(dx));
    }
  };

  template <typename RealElement, typename Basis, typename FERequirementTypeImpl = FErequirements<Eigen::VectorXd>>
  class AutoDiffFEClean {
  public:
    using LocalView = typename Basis::LocalView;
    using Traits    = TraitsFromLocalView<LocalView>;
    explicit AutoDiffFEClean(const Basis& basis, const typename LocalView::Element& element) {
      auto localView = basis.localView();
      localView.bind(element);
      localdofSize = localView.size();
    }
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      Eigen::VectorXdual2nd dx(localdofSize);
      Eigen::VectorXd g;
      autodiff::dual2nd e;
      dx.setZero();
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      hessian(f, wrt(dx), at(dx), e, g, h);
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      Eigen::VectorXdual dx(localdofSize);
      dx.setZero();
      autodiff::dual e;
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      gradient(f, wrt(dx), at(dx), e, g);
    }

    void calculateLocalSystem(const FERequirementType& par, typename Traits::MatrixType& h,
                              typename Traits::VectorType& g) const {
      Eigen::VectorXdual2nd dx(localdofSize);
      dx.setZero();
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      hessian(f, wrt(dx), at(dx), g, h);
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(localdofSize);
      dx.setZero();

      return this->underlying().calculateScalarImpl(par, dx);
    }

    size_t size() const { return localdofSize; }

  private:
    RealElement const& underlying() const  // CRTP
    {
      return static_cast<RealElement const&>(*this);
    }

    int localdofSize{};
  };

}  // namespace Ikarus
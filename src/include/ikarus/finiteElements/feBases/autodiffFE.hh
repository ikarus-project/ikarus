////
//// Created by Alex on 05.07.2021.
////
//
#pragma once
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include "../interface/finiteElementFunctionConcepts.hh"
#include "../physicsHelper.hh"

namespace Ikarus {
  ////This element can not be used on its own but it should be inherited from
  //// The class constructor can only be called from the templated class.
  template <typename RealElement, typename Basis, typename FERequirementTypeImpl = FErequirements<Eigen::VectorXd>>
  class AutoDiffFEClean {
  public:
    using LocalView = typename Basis::LocalView;
    using Traits    = TraitsFromLocalView<LocalView>;

    using FERequirementType = FErequirements<Eigen::VectorXd>;
    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      Eigen::VectorXdual2nd dx(localdofSize);
      Eigen::VectorXd g;
      autodiff::dual2nd e;
      dx.setZero();
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      Eigen::VectorXdual dx(localdofSize);
      dx.setZero();
      autodiff::dual e;
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      gradient(f, autodiff::wrt(dx), at(dx), e, g);
    }

    void calculateLocalSystem(const FERequirementType& par, typename Traits::MatrixType& h,
                              typename Traits::VectorType& g) const {
      Eigen::VectorXdual2nd dx(localdofSize);
      dx.setZero();
      auto f = [&](auto& x) { return this->underlying().calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), g, h);
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(localdofSize);
      dx.setZero();

      return this->underlying().calculateScalarImpl(par, dx);
    }

    size_t size() const { return localdofSize; }

  protected:
    explicit AutoDiffFEClean(const Basis& basis, const typename LocalView::Element& element) {
      auto localView = basis.localView();
      localView.bind(element);
      localdofSize = localView.size();
    }

  private:
    RealElement const& underlying() const  // CRTP
    {
      return static_cast<RealElement const&>(*this);
    }

    int localdofSize{};
  };

}  // namespace Ikarus
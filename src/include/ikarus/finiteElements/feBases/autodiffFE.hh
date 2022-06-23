/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

////
//
#pragma once
#include "../physicsHelper.hh"

#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feRequirements.hh>

namespace Ikarus {
  ////This element can not be used on its own but it should be inherited from
  //// The class constructor can only be called from the templated class.
  template <typename RealElement, typename Basis, typename FERequirementTypeImpl = FErequirements<Eigen::VectorXd>>
  class AutoDiffFE {
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
    explicit AutoDiffFE(const Basis& basis, const typename LocalView::Element& element) {
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
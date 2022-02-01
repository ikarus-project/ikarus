
// /*
//  *  This file is part of the Ikarus distribution (https://github.com/rath3t/Ikarus).
//  *  Copyright (c) 2021 Alexander Müller.
//  *  Institut fuer Baustatik und Baudynamik
//  *  Universität Stuttgart
//  *
//  *  This library is free software; you can redistribute it and/or
//  *   modify it under the terms of the GNU Lesser General Public
//  *   License as published by the Free Software Foundation; either
//  *   version 2.1 of the License, or (at your option) any later version.
//
// *   This library is distributed in the hope that it will be useful,
// *   but WITHOUT ANY WARRANTY; without even the implied warranty of
// *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// *   Lesser General Public License for more details.
//
// *   You should have received a copy of the GNU Lesser General Public
// *   License along with this library; if not, write to the Free Software
// *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// *  USA
// *

#pragma once
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <concepts>
#include <iostream>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <spdlog/spdlog.h>

#include "ikarus/AnsatzFunctions/Lagrange.h"
#include "ikarus/utils/LinearAlgebraHelper.h"
#include <ikarus/FiniteElements/FEPolicies.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Geometries/GeometryWithExternalInput.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>
#include <ikarus/Geometries/GeometryType.h>

namespace Ikarus::Variable {
  class IVariable;
}

namespace Ikarus::FiniteElements {

  template <typename LocalView>
  class NonLinearElasticityFEWithLocalBasis  {
  public:
//    using Base = FEVertexDisplacement<GridElementEntityType, IndexSetType>;
    NonLinearElasticityFEWithLocalBasis(LocalView& localView, double emod, double nu)
        :  localView_{localView}, emod_{emod}, nu_{nu} {}

    struct Traits{
      using GridEntity = typename LocalView::Element;
      /** \brief Dimension of the world space */
      static constexpr int worlddim = GridEntity::Geometry::coorddimension;

      /** \brief Dimension of the geometry */
      static constexpr int mydim = GridEntity::mydimension;

      /** \brief Dimension of the grid */
      static constexpr int dimension = GridEntity::dimension;

      /** \brief Type of the internal forces */
      using VectorType = Eigen::VectorXd;

      /** \brief Type of the stiffness matrix */
      using MatrixType = Eigen::MatrixXd;

      /** \brief Type of the FE parameters */
      using FERequirementType = typename IFiniteElement::FERequirementType;
    };

    using DeformedGeometry = Ikarus::Geometry::GeometryWithExternalInput<double, Traits::mydim, Traits::dimension>;

    [[nodiscard]] typename Traits::MatrixType calculateMatrix(const Eigen::VectorXd& displacements, const double& lambda) const {
      Eigen::VectorXdual2nd dx(localView_.size());
      dx.setZero();
      auto f = [&](Eigen::VectorXdual2nd &x) { return calculateScalarImpl<autodiff::dual2nd>(displacements,lambda,x); };
      const auto h = hessian(f, wrt(dx), at(dx));
      auto first_child = localView_.tree().child(0);
      Eigen::MatrixXd hE;
      hE.setZero(localView_.size(),localView_.size());
      for (auto i = 0U; i < localView_.size(); ++i) {
      for (auto j = 0U; j < localView_.size(); ++j)
        hE(i,j) = h(i,j);
      }
      return hE;
    }

    [[nodiscard]] typename Traits::VectorType calculateVector(const Eigen::VectorXd& displacements, const double& lambda) const {
      Eigen::VectorXdual dx(localView_.size());
      dx.setZero();
      auto f = [&](Eigen::VectorXdual &x) { return calculateScalarImpl<autodiff::dual>(displacements,lambda, x); };
      const auto g = gradient(f, wrt(dx), at(dx));
      Eigen::VectorXd gE;
      gE.setZero(localView_.size());
      for (auto i = 0U; i < localView_.size(); ++i) {
        gE[i] = g[i];
      }
      return gE;
    }

    template <class ScalarType>
    ScalarType calculateScalarImpl([[maybe_unused]] const Eigen::VectorXd& displacements, const double& lambda, Eigen::VectorX<ScalarType>& dx) const {
      Eigen::VectorX<ScalarType> localDisp(localView_.size());

      auto first_child = localView_.tree().child(0);
      const auto& fe = first_child.finiteElement();

        for(auto i = 0U; i< fe.size(); ++i)
          for(auto k2 = 0U; k2< Traits::mydim; ++k2)
        localDisp[2*i+k2] =  dx[localView_.tree().child(k2).localIndex(i)]+ displacements[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];
      ScalarType energy = 0.0;

      Eigen::Matrix<ScalarType, Traits::dimension, 4> dxM;
      dxM.setZero();
      for (auto pos = 0U, i = 0U; i < fe.size(); ++i) {
        dxM.col(i) = localDisp.template segment<Traits::mydim>(pos);
        pos += Traits::dimension;
      }
      const int order = 2*(fe.localBasis().order());
      const auto rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
      Eigen::Matrix3<ScalarType> C;
      C.setZero();
      const double fac = emod_ / (1 - nu_ * nu_);
      C(0, 0) = C(1, 1) = 1;
      C(0, 1) = C(1, 0) = nu_;
      C(2, 2)           = (1 - nu_) / 2;
      C *= fac;
      for (auto &gp : rule) {
        const auto geo = localView_.element().geometry();
        const auto J   = toEigenMatrix(geo.jacobianTransposed(gp.position()));
        //        std::cout<<"J:\n"<<J<<std::endl;
        Eigen::Matrix<ScalarType, Traits::dimension, 4> x;
        x.setZero();
        for (auto i = 0U; i < fe.size(); ++i)
          x.col(i) = dxM.col(i) + toEigenVector(geo.corner(i)) ;

        const Ikarus::Geometry::GeometryWithExternalInput<ScalarType, Traits::mydim, Traits::dimension> deformedgeo;

        std::vector<Dune::FieldMatrix<double,1,Traits::mydim>> dNM;
        fe.localBasis().evaluateJacobian(gp.position(),dNM);
        std::vector<Dune::FieldVector<double,1>> NM;
        fe.localBasis().evaluateFunction(gp.position(),NM);
        Eigen::Vector<ScalarType,Traits::worlddim> xv;
        xv.setZero();
        for (size_t i = 0; i < NM.size(); ++i) {
          xv+= x.col(i)*NM[i];
        }
        Eigen::Matrix<double, Eigen::Dynamic, Traits::mydim> dN;
        dN.resize(fe.size(),Eigen::NoChange);
        dN.setZero();
        for (auto i = 0U; i < fe.size(); ++i) {
          for (int j = 0; j < Traits::mydim; ++j) {
            dN(i,j) = dNM[i][0][j];
          }
        }
        const auto Jloc = Ikarus::LinearAlgebra::orthonormalizeMatrixColumns(J.transpose());
        const auto j    = deformedgeo.jacobianTransposed(dN, x);
        const Eigen::Matrix<ScalarType, Traits::mydim, Traits::mydim> F = Jloc.transpose() * j * (J.inverse()) * Jloc;
        dN *= (J * Jloc).inverse();

        const auto E
            = (0.5 * (F.transpose() * F - Eigen::Matrix<ScalarType, Traits::mydim, Traits::mydim>::Identity())).eval();
        Eigen::Vector<ScalarType, (Traits::mydim * (Traits::mydim + 1)) / 2> EVoigt;
        EVoigt.setZero();
        for (int i = 0; i < Traits::mydim; ++i)
          EVoigt(i) = E(i, i);

        if constexpr (Traits::mydim > 1) EVoigt(Traits::mydim) = E(0, 1) * 2;
        if constexpr (Traits::mydim > 2) {
          EVoigt(Traits::mydim + 1) = E(0, 2) * 2;
          EVoigt(Traits::mydim + 2) = E(1, 2) * 2;
        }
        Eigen::Vector<double,Traits::worlddim> fext;
        fext.setZero();
        fext[0] = 0;
        fext[1] = lambda;
        energy += EVoigt.dot(C * EVoigt) * geo.integrationElement(gp.position()) * gp.weight();
        energy -= xv.dot(fext) * geo.integrationElement(gp.position()) * gp.weight();
      }
      return energy;
    }

    LocalView localView_;
    double emod_;
    double nu_;
  };

}  // namespace Ikarus::FiniteElements

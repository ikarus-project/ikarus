
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

namespace Ikarus::Variable {
  class IVariable;
}

namespace Ikarus::FiniteElements {

  template <typename GridElementEntityType, typename IndexSetType>
  class NonLinearElasticityFE : public FEVertexDisplacement<GridElementEntityType, IndexSetType> {
  public:
    using Base = FEVertexDisplacement<GridElementEntityType, IndexSetType>;
    NonLinearElasticityFE(GridElementEntityType &gE, const IndexSetType &indexSet, double emod, double nu)
        : Base(gE, indexSet), elementGridEntity{&gE}, indexSet_{&indexSet}, emod_{emod}, nu_{nu} {}

    using Traits = FETraits<GridElementEntityType>;

    using DeformedGeometry = Ikarus::Geometry::GeometryWithExternalInput<double, Traits::mydim, Traits::dimension>;

    void initialize() {}

    [[nodiscard]] std::pair<typename Traits::MatrixType, typename Traits::VectorType> calculateLocalSystem(
        const typename Traits::FERequirementType &par) const {
      if (par.matrixAffordances == stiffness && par.vectorAffordances == forces)
        return calculateStiffnessMatrixAndInternalForcesImpl(par);
      else
        throw std::logic_error("This element can not handle your affordance! ");
    }

    [[nodiscard]] typename Traits::MatrixType calculateMatrix(const typename Traits::FERequirementType &req) const {
      if (req.matrixAffordances == stiffness)
        return calculateStiffnessMatrixAndInternalForcesImpl<false, true>(req);
      else
        throw std::logic_error("This element can not handle your affordance! ");
    }

    [[nodiscard]] double calculateScalar([[maybe_unused]] const typename Traits::FERequirementType &req) const {
      if (req.scalarAffordances == potentialEnergy) {
        Eigen::VectorXd d(this->dofSize());
        auto dv = req.variables.value().get().get(EntityType::vertex);
        for (int pos = 0, i = 0; i < this->dofSize() / 2; ++i) {
          d.template segment<Traits::mydim>(pos) = Variable::getValue(dv[i]);
          pos += Traits::mydim;
        }
        return 0.5 * d.dot(calculateStiffnessMatrixAndInternalForcesImpl<false, true>(req) * d);
      } else
        throw std::logic_error("This element can not handle your scalar affordance! ");
    }

    [[nodiscard]] typename Traits::VectorType calculateVector(const typename Traits::FERequirementType &req) const {
      return calculateStiffnessMatrixAndInternalForcesImpl<true, false>(req);
    }

    template <bool internalForcesFlag = true, bool stiffnessMatrixFlag = true>
    auto calculateStiffnessMatrixAndInternalForcesImpl(
        [[maybe_unused]] const typename Traits::FERequirementType &req) const {
      if constexpr (internalForcesFlag && stiffnessMatrixFlag) {
        const auto rule = Dune::QuadratureRules<double, Traits::mydim>::rule(duneType(elementGridEntity->type()), 2);
        typename Traits::VectorType Fint(this->dofSize());
        typename Traits::MatrixType K(this->dofSize(), this->dofSize());
        Fint.setZero();
        K.setZero();
        return std::make_pair(K, Fint);
      } else if constexpr (internalForcesFlag && !stiffnessMatrixFlag) {
        typename Traits::VectorType Fint(this->dofSize());
        Fint.setZero();
        const auto rule = Dune::QuadratureRules<double, Traits::mydim>::rule(duneType(elementGridEntity->type()), 2);
        Eigen::Matrix3d C;
        C.setZero();
        const double fac = emod_ / (1 - nu_ * nu_);
        C(0, 0) = C(1, 1) = 1;
        C(0, 1) = C(1, 0) = nu_;
        C(2, 2)           = (1 - nu_) / 2;
        C *= fac;
        for (auto &gp : rule) {
          const auto geo = elementGridEntity->geometry();
          const auto J   = toEigenMatrix(geo.jacobianTransposed(gp.position()));
          Eigen::Matrix<double, Traits::dimension, 4> x;

          auto dv = req.variables.value().get().get(EntityType::vertex);
          for (int pos = 0, i = 0; i < 4; ++i) {
            x.col(i) = Variable::getValue(dv[i]) + toEigenVector(geo.corner(i));
          }
          const DeformedGeometry deformedgeo;

          Eigen::Matrix<double, 4, Traits::mydim> dN
              = Ikarus::LagrangeCube<double, Traits::mydim, 1>::evaluateJacobian(toEigenVector(gp.position()));
          const auto Jloc = Ikarus::LinearAlgebra::orthonormalizeMatrixColumns(J.transpose());
          //          dN *= (J * Jloc).inverse();
          const auto j                                                = deformedgeo.jacobianTransposed(dN, x);
          const Eigen::Matrix<double, Traits::mydim, Traits::mydim> F = Jloc.transpose() * j * (J.inverse()) * Jloc;
                    dN *= (J * Jloc).inverse();
          const auto bop = boperator(dN, F);
//          std::cout<<"FinFint:\n"<<F<<std::endl;
          const auto E
              = (0.5 * (F.transpose() * F - Eigen::Matrix<double, Traits::mydim, Traits::mydim>::Identity())).eval();
          Eigen::Vector<double, (Traits::mydim * (Traits::mydim + 1)) / 2> EVoigt;
          for (int i = 0; i < Traits::mydim; ++i)
            EVoigt(i) = E(i, i);

          if constexpr (Traits::mydim > 1) EVoigt(Traits::mydim) = E(0, 1) * 2;
          if constexpr (Traits::mydim > 2) {
            EVoigt(Traits::mydim + 1) = E(0, 2) * 2;
            EVoigt(Traits::mydim + 2) = E(1, 2) * 2;

          }
//          std::cout<<"EinFint:\n"<<E<<std::endl;
//          std::cout<<"EVoigtinFint:\n"<<EVoigt<<std::endl;
          Fint += (bop.transpose() * C * EVoigt) * geo.integrationElement(gp.position()) * gp.weight();
        }
//        std::cout << "F:\n" << Fint.transpose() << std::endl;
        return Fint;
      } else if constexpr (not internalForcesFlag && stiffnessMatrixFlag) {
        typename Traits::MatrixType K(this->dofSize(), this->dofSize());
        K.setZero();
        const auto rule = Dune::QuadratureRules<double, Traits::mydim>::rule(duneType(elementGridEntity->type()), 2);
        Eigen::Matrix3d C;
        C.setZero();
        const double fac = emod_ / (1 - nu_ * nu_);
        C(0, 0) = C(1, 1) = 1;
        C(0, 1) = C(1, 0) = nu_;
        C(2, 2)           = (1 - nu_) / 2;
        C *= fac;
        for (auto &gp : rule) {
          const auto geo = elementGridEntity->geometry();
          const auto J   = toEigenMatrix(geo.jacobianTransposed(gp.position()));
          Eigen::Matrix<double, Traits::dimension, 4> x;

          auto dv = req.variables.value().get().get(EntityType::vertex);
          for (int i = 0; i < 4; ++i) {
            x.col(i) = Variable::getValue(dv[i]) + toEigenVector(geo.corner(i));
          }
          std::cout<<"x:\n"<<x<<std::endl;
          const DeformedGeometry deformedgeo;

          Eigen::Matrix<double, 4, Traits::mydim> dN
              = Ikarus::LagrangeCube<double, Traits::mydim, 1>::evaluateJacobian(toEigenVector(gp.position()));
          const auto Jloc = Ikarus::LinearAlgebra::orthonormalizeMatrixColumns(J.transpose());
          //          dN *= (J * Jloc).inverse();
          const auto j                                                = deformedgeo.jacobianTransposed(dN, x);
          const Eigen::Matrix<double, Traits::mydim, Traits::mydim> F = Jloc.transpose() * j * (J.inverse()) * Jloc;
                    dN *= (J * Jloc).inverse();
          const auto bop = boperator(dN, F);
//                    std::cout<<"F:\n"<<F<<std::endl;
          //          std::cout<<"B:\n"<<bop<<std::endl;
          //          std::cout<<"C:\n"<<C<<std::endl;
          K += (bop.transpose() * C * bop) * geo.integrationElement(gp.position()) * gp.weight();
        }
        //        std::cout<<"K:\n"<<K<<std::endl;
        return K;
      } else
        static_assert(internalForcesFlag == false && stiffnessMatrixFlag == false,
                      "You asked the element: \"Don't return anything\"");
    }

    [[nodiscard]] unsigned int subEntities(unsigned int codim) const { return elementGridEntity->subEntities(codim); }
    [[nodiscard]] unsigned int subIndex(int i, unsigned int codim) const {
      return indexSet_->subIndex(*elementGridEntity, i, codim);
    }

    [[nodiscard]] unsigned int dimension() const { return Traits::mydim; }

  private:
    auto boperator(const Eigen::Matrix<double, 4, Traits::mydim> &dN,
                   const Eigen::Matrix<double, Traits::mydim, Traits::mydim> &F) const {
      typename Traits::MatrixType Bop((Traits::mydim * (Traits::mydim + 1)) / 2, this->dofSize());
      Bop.setZero();
      for (int i = 0; i < dN.rows(); ++i) {
        auto currentBlock
            = Bop.template block<(Traits::mydim * (Traits::mydim + 1)) / 2, Traits::mydim>(0, Traits::mydim * i);

        for (int j = 0; j < Traits::mydim; ++j)
          currentBlock.template block<1, Traits::mydim>(j, 0) = (F.col(j) * dN(i, j)).transpose();

        if constexpr (Traits::mydim > 1)
          currentBlock.template block<1, Traits::mydim>(Traits::mydim, 0)
              = (F.col(0) * dN(i, 1) + F.col(1) * dN(i, 0)).transpose();
        if constexpr (Traits::mydim > 2) {
          currentBlock.template block<1, Traits::mydim>(4, 0) = (F.col(1) * dN(i, 2) + F.col(2) * dN(i, 1)).transpose();
          currentBlock.template block<1, Traits::mydim>(5, 0) = (F.col(0) * dN(i, 2) + F.col(2) * dN(i, 0)).transpose();
        }
      }
      return Bop;
    }
    GridElementEntityType const *const elementGridEntity;
    IndexSetType const *const indexSet_;
    double emod_;
    double nu_;
  };

}  // namespace Ikarus::FiniteElements

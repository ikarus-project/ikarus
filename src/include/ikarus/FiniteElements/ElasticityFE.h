
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
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {
  class IVariable;
}

namespace Ikarus::FiniteElements {

  template <typename GridElementEntityType, typename IndexSetType>
  class ElasticityFE : public FEVertexDisplacement<GridElementEntityType, IndexSetType> {
  public:
    using Base = FEVertexDisplacement<GridElementEntityType, IndexSetType>;
    ElasticityFE(GridElementEntityType &gE, const IndexSetType &indexSet, double emod, double nu)
        : Base(gE, indexSet), elementGridEntity{gE}, indexSet_{&indexSet}, emod_{emod}, nu_{nu} {}

    using Traits = FETraits<GridElementEntityType>;

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
        const auto rule = Dune::QuadratureRules<double, Traits::mydim>::rule(duneType(elementGridEntity.type()), 2);
        typename Traits::VectorType Fint(this->dofSize());
        typename Traits::MatrixType K(this->dofSize(), this->dofSize());
        Fint.setZero();
        K.setZero();
        return std::make_pair(K, Fint);
      } else if constexpr (internalForcesFlag && !stiffnessMatrixFlag) {
        Eigen::VectorXd d(this->dofSize());
        auto dv = req.variables.value().get().get(EntityType::vertex);
        for (int pos = 0, i = 0; i < this->dofSize() / 2; ++i) {
          d.template segment<Traits::mydim>(pos) = Variable::getValue(dv[i]);
          pos += Traits::mydim;
        }
        Eigen::VectorXd fint = calculateStiffnessMatrixAndInternalForcesImpl<false, true>(req) * d;
        return fint;
      } else if constexpr (not internalForcesFlag && stiffnessMatrixFlag) {
        typename Traits::MatrixType K(this->dofSize(), this->dofSize());
        K.setZero();
        const auto rule = Dune::QuadratureRules<double, Traits::mydim>::rule(duneType(elementGridEntity.type()), 2);
        Eigen::Matrix3d C;
        C.setZero();
        const double fac = emod_ / (1 - nu_ * nu_);
        C(0, 0) = C(1, 1) = 1;
        C(0, 1) = C(1, 0) = nu_;
        C(2, 2)           = (1 - nu_) / 2;
        C *= fac;
        for (auto &gp : rule) {
          const auto geo = elementGridEntity.geometry();
          const auto J   = toEigenMatrix(geo.jacobianTransposed(gp.position()));
          Eigen::Matrix<double, 4, Traits::mydim> dN
              = Ikarus::LagrangeCube<double, Traits::mydim, 1>::evaluateJacobian(toEigenVector(gp.position()));
          dN *= (J * Ikarus::LinearAlgebra::orthonormalizeMatrixColumns(J.transpose())).inverse();
          const auto bop = boperator(dN);
          K += (bop.transpose() * C * bop) * geo.integrationElement(gp.position()) * gp.weight();
        }
        return K;
      } else
        static_assert(internalForcesFlag == false && stiffnessMatrixFlag == false,
                      "You asked the element: \"Don't return anything\"");
    }

    [[nodiscard]] unsigned int subEntities(unsigned int codim) const { return elementGridEntity.subEntities(codim); }
    [[nodiscard]] unsigned int subIndex(int i, unsigned int codim) const {
      return indexSet_->subIndex(elementGridEntity, i, codim);
    }

    [[nodiscard]] unsigned int dimension() const { return Traits::mydim; }

  private:
    auto boperator(const Eigen::Matrix<double, 4, Traits::mydim> &dN) const {
      typename Traits::MatrixType Bop((Traits::mydim * (Traits::mydim + 1)) / 2, this->dofSize());
      Bop.setZero();
      for (int i = 0; i < dN.rows(); ++i) {
        for (int j = 0; j < Traits::mydim; ++j)
          Bop(j, Traits::mydim * i + j) = dN(i, j);

        for (int s = 0; s < Traits::mydim * (Traits::mydim - 1) / 2; ++s)  // loop over off-diagonal strains
        {
          std::array<int, 2> curdim{};
          for (int c = 0, d = 0; d < Traits::mydim; ++d) {
            if (d == 2 - s) continue;
            curdim[c++] = d;
          }
          for (int k = 0; k < 2; ++k)
            Bop(s + Traits::mydim, Traits::mydim * i + curdim[k]) = dN(i, curdim[(k + 1) % 2]);
          //              Bop(s + mydimension, mydimension * i + curdim[1]) = dN(i, curdim[0]);
        }
      }
      return Bop;
    }
    GridElementEntityType elementGridEntity;
    IndexSetType const *const indexSet_;
    double emod_;
    double nu_;
  };

}  // namespace Ikarus::FiniteElements

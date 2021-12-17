
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
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {
  class IVariable;
}

namespace Ikarus::FiniteElements {

  template <typename GridElementEntityType, typename IndexSetType, std::floating_point ct = double>
  class ElasticityFE {
  public:
    ElasticityFE(GridElementEntityType &gE, const IndexSetType &indexSet, double emod, double nu)
        : elementGridEntity{&gE}, indexSet_{&indexSet}, emod_{emod}, nu_{nu} {}

    /** \brief Type used for coordinates */
    using ctype = ct;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = GridElementEntityType::dimensionworld;

    /** \brief Dimension of the geometry */
    static constexpr int mydim = GridElementEntityType::mydimension;

    /** \brief Type of the Nodes coordinate */
    using NodeType = Eigen::Matrix<ctype, worlddim, 1>;

    /** \brief Type of the ParameterSpace coordinate */
    using ParameterSpaceType = Eigen::Matrix<ctype, mydim, 1>;

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using DofTupleVectorType = typename IFiniteElement::DofPairVectorType;

    /** \brief Type of the FE parameters */
    using FERequirementType = typename IFiniteElement::FERequirementType;

    /** \brief Type of the Variables */
    using VariableVectorType = typename FERequirementType::VariableType;

    /** \brief Type of the DataVector */
    using DataVectorType = typename FERequirementType::DataType ;

    /** \brief Type of the Dofs / SolutionType
     * using NodalSolutionType = Displacement<ctype,worlddim>;*/

    /** \brief Type of the internal forces */
    using VectorType = Eigen::VectorXd;

    /** \brief Type of the stiffness matrix */
    using MatrixType = Eigen::MatrixXd;

    [[nodiscard]] constexpr int dofSize() const {
      if constexpr (worlddim == 3)
        return vertices(*elementGridEntity).size() * 3;
      else if constexpr (worlddim == 2)
        return vertices(*elementGridEntity).size() * 2;
    }

    void initialize() {}

    [[nodiscard]] std::pair<MatrixType, VectorType> calculateLocalSystem(const FERequirementType & par) const {
      if (par.matrixAffordances == stiffness && par.vectorAffordances == forces)
        return calculateStiffnessMatrixAndInternalForcesImpl(par);
      else
        throw std::logic_error("This element can not handle your affordance! ");
    }

    [[nodiscard]] MatrixType calculateMatrix(const FERequirementType & par) const {
      if (par.matrixAffordances == stiffness)
        return calculateStiffnessMatrixAndInternalForcesImpl<false, true>(par);
      else
        throw std::logic_error("This element can not handle your affordance! ");
    }

    [[nodiscard]] double calculateScalar([[maybe_unused]]  const FERequirementType & par) const {
      return 13.0;
    }

    [[nodiscard]] VectorType calculateVector(const FERequirementType & par) const {
      return calculateStiffnessMatrixAndInternalForcesImpl<true, false>(par);
    }

    template <bool internalForcesFlag = true, bool stiffnessMatrixFlag = true>
    auto calculateStiffnessMatrixAndInternalForcesImpl([[maybe_unused]] const FERequirementType &req) const {

      if constexpr (internalForcesFlag && stiffnessMatrixFlag) {
        const auto rule = Dune::QuadratureRules<ctype, mydim>::rule(duneType(elementGridEntity->type()), 2);

        VectorType Fint(dofSize());
        MatrixType K(dofSize(), dofSize());
        Fint.setZero();
        K.setZero();
        for (auto &gp : rule) {
          const auto dN         = Ikarus::LagrangeCube<double, 2, 1>::evaluateJacobian(toEigenVector(gp.position()));
          const auto geo        = elementGridEntity->geometry();
          static_assert(std::remove_cvref_t<decltype(gp.position())>::dimension == 2);
          const auto refJacobian = geo.jacobianInverseTransposed(gp.position());
          //          Eigen::Vector2d u
          //          vol += integrationElement(gp.position()) * gp.weight();
        }

        return std::make_pair(K, Fint);
      } else if constexpr (internalForcesFlag && !stiffnessMatrixFlag) {
        return VectorType::Ones(8);
      } else if constexpr (not internalForcesFlag && stiffnessMatrixFlag) {
        VectorType Fint(dofSize());
        MatrixType K(dofSize(), dofSize());
        MatrixType Bop((mydim * (mydim + 1)) / 2, dofSize());
        K.setZero();
        Bop.setZero();
        const auto rule = Dune::QuadratureRules<ctype, mydim>::rule(duneType(elementGridEntity->type()), 2);
        double vol      = 0;
        Eigen::Matrix3d C;
        C.setZero();
        const double fac = emod_ / (1 - nu_ * nu_);
        C(0, 0) = C(1, 1) = 1;
        C(0, 1) = C(1, 0) = nu_;
        C(2, 2)           = (1 - nu_) / 2;
        C *= fac;
        for (auto &gp : rule) {
          const auto geo        = elementGridEntity->geometry();
          const auto J          = toEigenMatrix(geo.jacobianTransposed(gp.position()));
          Eigen::Matrix<double, 4, mydim> dN
              = Ikarus::LagrangeCube<double, mydim, 1>::evaluateJacobian(toEigenVector(gp.position()));
          dN *= (J * Ikarus::LinearAlgebra::orthonormalizeMatrixColumns(J.transpose())).inverse();

          for (int i = 0; i < dN.rows(); ++i) {
            for (int j = 0; j < mydim; ++j)
              Bop(j, mydim * i + j) = dN(i, j);

            for (int s = 0; s < mydim * (mydim - 1) / 2; ++s)  // loop over off-diagonal strains
            {
              std::array<int, 2> curdim{};
              for (int c = 0, d = 0; d < mydim; ++d) {
                if (d == 2 - s) continue;
                curdim[c++] = d;
              }
              for (int k = 0; k < 2; ++k)
                Bop(s + mydim, mydim * i + curdim[k]) = dN(i, curdim[(k + 1) % 2]);
              //              Bop(s + mydimension, mydimension * i + curdim[1]) = dN(i, curdim[0]);
            }
          }
          K += (Bop.transpose() * C * Bop) * geo.integrationElement(gp.position()) * gp.weight();
          vol += geo.integrationElement(gp.position()) * gp.weight();
        }
        return K;
      } else
        static_assert(internalForcesFlag == false && stiffnessMatrixFlag == false,
                      "You asked the element: \"Don't return anything\"");
    }

    [[nodiscard]] DofTupleVectorType getEntityVariableTuple() const {
      DofTupleVectorType entDofTupleVector(vertices(*elementGridEntity).size());
      using namespace Ikarus::Variable;
      VariableTags dofType;
      if constexpr (worlddim == 3)
        dofType = VariableTags::displacement3d;
      else if constexpr (worlddim == 2)
        dofType = VariableTags::displacement2d;
      else if constexpr (worlddim == 1)
        dofType = VariableTags::displacement1d;
      else
        static_assert(worlddim > 3 || worlddim < 1, "This element has an impossible worlddim.");
      for (int id = 0; auto &entityDofTuple : entDofTupleVector) {
        entityDofTuple.entityID = indexSet_->subIndex(*elementGridEntity, id++, worlddim);
        entityDofTuple.variableVector.assign(1, dofType);
        entityDofTuple.entityType = EntityType::vertex;
      }

      return entDofTupleVector;
    }

    [[nodiscard]] unsigned int subEntities(unsigned int codim) const { return elementGridEntity->subEntities(codim); }
    [[nodiscard]] unsigned int subIndex(int i, unsigned int codim) const {
      return indexSet_->subIndex(*elementGridEntity, i, codim);
    }

    [[nodiscard]] unsigned int dimension() const { return mydim; }

  private:
    GridElementEntityType const *const elementGridEntity;
    IndexSetType const *const indexSet_;
    double emod_;
    double nu_;
  };

}  // namespace Ikarus::FiniteElements


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
#include <dune/geometry/type.hh>

#include <spdlog/spdlog.h>

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
    ElasticityFE(GridElementEntityType &gE, const IndexSetType &indexSet)
        : elementGridEntity{&gE}, indexSet_{&indexSet} {}

    /** \brief Type used for coordinates */
    using ctype = ct;

    /** \brief Dimension of the world space */
    static constexpr int coorddimension = GridElementEntityType::dimensionworld;

    /** \brief Dimension of the geometry */
    static constexpr int mydimension = GridElementEntityType::mydimension;

    /** \brief Type of the Nodes coordinate */
    using NodeType = Eigen::Matrix<ctype, coorddimension, 1>;

    /** \brief Type of the ParameterSpace coordinate */
    using ParameterSpaceType = Eigen::Matrix<ctype, mydimension, 1>;

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using DofTupleVectorType = typename IFiniteElement::DofPairVectorType;

    /** \brief Type of the DofVector */
    using VariableVectorType = typename IFiniteElement::VariableVectorType;

    /** \brief Type of the DataVector */
    using DataVectorType = typename IFiniteElement::DataVectorType;

    /** \brief Type of the Dofs / SolutionType
     * using NodalSolutionType = Displacement<ctype,coorddimension>;*/

    /** \brief Type of the internal forces */
    using VectorType = Eigen::VectorXd;

    /** \brief Type of the stiffness matrix */
    using MatrixType = Eigen::MatrixXd;

    [[nodiscard]] constexpr int dofSize() const {
      if constexpr (coorddimension == 3)
        return vertices(*elementGridEntity).size() * 3;
      else if constexpr (coorddimension == 2)
        return vertices(*elementGridEntity).size() * 2;
    }

    void initialize() {}

    [[nodiscard]] std::pair<MatrixType, VectorType> calculateLocalSystem(const MatrixAffordances &matA,
                                                                         const VectorAffordances &vecA,
                                                                         [[maybe_unused]] VariableVectorType &vars,
                                                                         [[maybe_unused]] DataVectorType data
                                                                         = std::nullopt) const {
      if (matA == stiffness && vecA == forces)
        return calculateStiffnessMatrixAndInternalForcesImpl();
      else
        throw std::logic_error("This element can not handle your affordance! ");
    }

    [[nodiscard]] MatrixType calculateMatrix(const MatrixAffordances &matA, [[maybe_unused]] VariableVectorType &vars,
                                             [[maybe_unused]] const DataVectorType &data = std::nullopt) const {
      if (matA == stiffness)
        return calculateStiffnessMatrixAndInternalForcesImpl<false, true>();
      else
        throw std::logic_error("This element can not handle your affordance! ");
    }

    [[nodiscard]] double calculateScalar(const ScalarAffordances &, [[maybe_unused]] VariableVectorType &vars,
                                         [[maybe_unused]] const DataVectorType &data = std::nullopt) const {
      return 13.0;
    }

    [[nodiscard]] VectorType calculateVector(const VectorAffordances &, [[maybe_unused]] VariableVectorType &vars,
                                             [[maybe_unused]] const DataVectorType &data = std::nullopt) const {
      return calculateStiffnessMatrixAndInternalForcesImpl<true, false>();
    }

    template <bool internalForcesFlag = true, bool stiffnessMatrixFlag = true>
    auto calculateStiffnessMatrixAndInternalForcesImpl() const {
      if constexpr (internalForcesFlag && stiffnessMatrixFlag) {
        const VectorType Fint = VectorType::Ones(8);
        const MatrixType K    = MatrixType::Ones(8, 8);
        return std::make_pair(K, Fint);
      } else if constexpr (internalForcesFlag && !stiffnessMatrixFlag)
        return VectorType::Ones(8);
      else if constexpr (!internalForcesFlag && stiffnessMatrixFlag)
        return MatrixType::Ones(8, 8);
      else
        static_assert(internalForcesFlag == false && stiffnessMatrixFlag == false,
                      "You asked the element: \"Don't return anything\"");
    }

    [[nodiscard]] DofTupleVectorType getEntityVariableTuple() const {
      DofTupleVectorType entDofTupleVector(vertices(*elementGridEntity).size());
      using namespace Ikarus::Variable;
      VariableTags dofType;
      if constexpr (coorddimension == 3)
        dofType = VariableTags::displacement3d;
      else if constexpr (coorddimension == 2)
        dofType = VariableTags::displacement2d;
      else if constexpr (coorddimension == 1)
        dofType = VariableTags::displacement1d;
      else
        static_assert(coorddimension > 3 || coorddimension < 1, "This element has an impossible coorddimension.");
      for (int id = 0; auto &entityDofTuple : entDofTupleVector) {
        entityDofTuple.entityID = indexSet_->subIndex(*elementGridEntity, id++, coorddimension);
        entityDofTuple.variableVector.assign(1, dofType);
        entityDofTuple.entityType = EntityType::vertex;
      }

      return entDofTupleVector;
    }

    [[nodiscard]] unsigned int subEntities(unsigned int codim) const { return elementGridEntity->subEntities(codim); }
    [[nodiscard]] unsigned int subIndex(int i, unsigned int codim) const {
      return indexSet_->subIndex(*elementGridEntity, i, codim);
    }

    [[nodiscard]] unsigned int dimension() const { return mydimension; }

  private:
    GridElementEntityType const *const elementGridEntity;
    IndexSetType const *const indexSet_;
  };

}  // namespace Ikarus::FiniteElements


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

#include <ikarus/FiniteElements/FiniteElementPolicies.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {
  class IVariable;
}

namespace Ikarus::FiniteElements {

  template <typename GridElementEntityType, std::floating_point ct = double>
  class ElasticityFE {
  public:
    explicit ElasticityFE(GridElementEntityType& gE) : elementGridEntity{&gE} {}

    /** \brief Type used for coordinates */
    using ctype = ct;

    /** \brief Dimension of the world space */
    static constexpr int coorddimension = GridElementEntityType::dimensionworld;

    /** \brief Dimension of the geometry */
    static constexpr int mydimension = GridElementEntityType::mydimension;

    /** \brief Type of the Nodes */
    using NodeType = Eigen::Matrix<ctype, coorddimension, 1>;

    /** \brief Type of the ParameterSpace coordinate */
    using ParameterSpaceType = Eigen::Matrix<ctype, mydimension, 1>;

    /** \brief Type of the DofVector */
    using DofVectorType = typename IFiniteElement::DofVectorType;

    /** \brief Type of the Dofs / SolutionType
     * using NodalSolutionType = Displacement<ctype,coorddimension>;*/

    /** \brief Type of the internal forces */
    using VectorType = DynVectord;

    /** \brief Type of the stiffness matrix */
    using MatrixType = DynMatrixd;

    [[nodiscard]] constexpr int dofSize() const {
      if constexpr (coorddimension == 3)
        return vertices(elementGridEntity).size() * 3;
      else if constexpr (coorddimension == 2)
        return vertices(elementGridEntity).size() * 2;
    }

    void initialize() { std::cout << "initialize ElasticityFE" << std::endl; }

    [[nodiscard]] std::pair<MatrixType, VectorType> calculateLocalSystem(const ElementMatrixAffordances& matA,
                                                                         const ElementVectorAffordances& vecA) const {
      if (matA == stiffness && vecA == forces)
        return calculateStiffnessMatrixAndInternalForcesImpl();
      else
        throw std::logic_error("This element can not handle your affordance! ");
    }

    [[nodiscard]] MatrixType calculateMatrix(const ElementMatrixAffordances& matA) const {
      if (matA == stiffness)
        return calculateStiffnessMatrixAndInternalForcesImpl<false, true>();
      else
        throw std::logic_error("This element can not handle your affordance! ");
    }

    [[nodiscard]] double calculateScalar(const ElementScalarAffordances&) const { return 13.0; }

    [[nodiscard]] VectorType calculateVector(const ElementVectorAffordances&) const {
      return calculateStiffnessMatrixAndInternalForcesImpl<true, false>();
    }

    template <bool internalForcesFlag = true, bool stiffnessMatrixFlag = true>
    auto calculateStiffnessMatrixAndInternalForcesImpl() const {
      if constexpr (internalForcesFlag && stiffnessMatrixFlag) {
        const VectorType Fint = VectorType::Ones(8);
        const MatrixType K    = MatrixType::Ones(8, 8);
        return std::make_pair(K, Fint);
      }

      else if constexpr (internalForcesFlag && !stiffnessMatrixFlag)
        return VectorType::Ones(8);
      else if constexpr (!internalForcesFlag && stiffnessMatrixFlag)
        return MatrixType::Ones(8, 8);
      else
        static_assert(internalForcesFlag == false && stiffnessMatrixFlag == false,
                      "You asked the element: \"Don't return anything\"");
    }

    [[nodiscard]] DofVectorType getEntityVariablePairs() const {
      DofVectorType dofs;
      for (auto&& vertex : vertices(elementGridEntity)) {
        if constexpr (coorddimension == 3)
          dofs.push_back({vertex->getID(), {Ikarus::Variable::displacement3d}});
        else if constexpr (coorddimension == 2)
          dofs.push_back({vertex->getID(), {Ikarus::Variable::displacement2d}});
      }

      return dofs;
    }

  private:
    GridElementEntityType const* const elementGridEntity;
  };

}  // namespace Ikarus::FiniteElements

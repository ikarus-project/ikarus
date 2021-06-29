
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

#include <spdlog/spdlog.h>

//#include <ikarus/Variables/GenericVariable.h>
#include <dune/geometry/type.hh>

#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {
  class GenericVariable;
}

namespace Ikarus::FiniteElements {
  template <typename GridElementEntityType, std::floating_point ct = double>
  class ElasticityFE {
  public:
    explicit ElasticityFE(GridElementEntityType& gE) : elementGridEntity{&gE} {
      if constexpr (coorddimension == 3)
        for (auto&& vertex : vertices(elementGridEntity))
          vertex->addDof(Ikarus::Variable::DISPLACEMENT3D());
    }

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
    using DofVectorType = std::vector<std::pair<size_t, std::vector<Ikarus::Variable::VariablesTags>>>;

    /** \brief Type of the Dofs / SolutionType
     * using NodalSolutionType = Displacement<ctype,coorddimension>;*/

    /** \brief Type of the internal forces */
    using VectorType = DynVectord;

    /** \brief Type of the stiffness matrix */
    using MatrixType = DynMatrixd;

    [[nodiscard]] constexpr int dofSize() const { return 3; }

    void generateDofs() const {
      //      std::vector<std::pair<
      //      for(auto&& vertex : vertices(elementGridEntity))
    }

    void initialize() { std::cout << "initialize ElasticityFE" << std::endl; }

    [[nodiscard]] std::pair<VectorType, MatrixType> calculateLocalSystem() const {
      return calculateStiffnessMatrixAndInternalForcesImpl();
    }

    [[nodiscard]] MatrixType calculateLHS() const {
      return calculateStiffnessMatrixAndInternalForcesImpl<false, true>();
    }

    [[nodiscard]] double getEnergy() const { return 0.0; }

    [[nodiscard]] VectorType calculateRHS() const {
      return calculateStiffnessMatrixAndInternalForcesImpl<true, false>();
    }

    template <bool internalForcesFlag = true, bool stiffnessMatrixFlag = true>
    auto calculateStiffnessMatrixAndInternalForcesImpl() const {
      if constexpr (internalForcesFlag && stiffnessMatrixFlag) {
        const VectorType Fint = VectorType::Ones(5);
        const MatrixType K    = MatrixType::Ones(5, 5);
        return std::make_pair(Fint, K);
      }

      else if constexpr (internalForcesFlag && !stiffnessMatrixFlag)
        return VectorType::Ones(5);
      else if constexpr (!internalForcesFlag && stiffnessMatrixFlag)
        return MatrixType::Ones(5, 5);
      else
        static_assert(internalForcesFlag == false && stiffnessMatrixFlag == false,
                      "You asked the element: \"Don't return anything\"");
    }

    [[nodiscard]] DofVectorType getDofVector() {
      DofVectorType dofs;

      for (auto&& vertex : vertices(elementGridEntity)) {
        if constexpr (coorddimension == 3)
          dofs.push_back({vertex->getID(), {Ikarus::Variable::displacement3d}});
        else if constexpr (coorddimension == 2)
          dofs.push_back({vertex->getID(), {Ikarus::Variable::displacement2d}});
      }

      return dofs;
    }

    auto getGridEntity() { return elementGridEntity; }

  private:
    GridElementEntityType const* const elementGridEntity;
  };

}  // namespace Ikarus::FiniteElements

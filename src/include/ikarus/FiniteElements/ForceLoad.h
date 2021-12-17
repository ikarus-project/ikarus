//
// Created by lex on 17/12/2021.
//

#pragma once
#include <concepts>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
namespace Ikarus::FiniteElements {

    template <typename GridElementEntityType, typename IndexSetType, std::floating_point ct = double>
    class ForceLoad {
    public:
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
      using DataVectorType = typename FERequirementType::DataType;

      /** \brief Type of the Dofs / SolutionType
     * using NodalSolutionType = Displacement<ctype,worlddim>;*/

      /** \brief Type of the internal forces */
      using VectorType = Eigen::VectorXd;

      /** \brief Type of the stiffness matrix */
      using MatrixType = Eigen::MatrixXd;

      ForceLoad(GridElementEntityType &gE, const IndexSetType &indexSet)
          : elementGridEntity{&gE}, indexSet_{&indexSet} {}

      [[nodiscard]] VectorType calculateVector(const FERequirementType &par) const { return calculateVectorImpl(par); }

      [[nodiscard]] VectorType calculateVectorImpl([[maybe_unused]] const FERequirementType &req) const {
        assert(req.parameter.contains(FEParameter::loadfactor));
        return VectorType{};
      }

    private:
      GridElementEntityType const *const elementGridEntity;
      IndexSetType const *const indexSet_;
    };


}
//
// Created by lex on 18/12/2021.
//

#pragma once
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>

namespace Ikarus::FiniteElements {

  template <typename GridElementEntityType>
  struct FETraits {
    /** \brief Type used for coordinates */
    using ctype = double;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = GridElementEntityType::Geometry::coorddimension;

    /** \brief Dimension of the geometry */
    static constexpr int mydim = GridElementEntityType::mydimension;

    /** \brief Dimension of the grid */
    static constexpr int dimension = GridElementEntityType::dimension;

    /** \brief Type of the  coordinate */
    using GlobalCoordinates = Eigen::Matrix<ctype, worlddim, 1>;

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
  };

  template <typename GridElementEntityType, typename IndexSetType>
  class FEVertexDisplacement {
  public:
    FEVertexDisplacement(GridElementEntityType &gE, const IndexSetType &indexSet)
        : elementGridEntity{gE}, indexSet_{&indexSet} {}

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using DofTupleVectorType = typename IFiniteElement::DofPairVectorType;

    using Traits = FETraits<GridElementEntityType>;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = Traits::worlddim;

    [[nodiscard]] constexpr int dofSize() const { return elementGridEntity.subEntities(Traits::dimension) * worlddim; }

    [[nodiscard]] DofTupleVectorType getEntityVariableTuple() const {
      DofTupleVectorType entDofTupleVector(elementGridEntity.subEntities(Traits::dimension));
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
        if constexpr (requires{ elementGridEntity->template subEntity<Traits::dimension>(id);})
          std::cout<<"id: "<<id<<" "<<"Coords: "<<elementGridEntity.template subEntity<Traits::dimension>(id).geometry().corner(0)<<" globID "<<indexSet_->subIndex(*elementGridEntity, id, Traits::dimension)<<std::endl;
        entityDofTuple.entityID = indexSet_->subIndex(elementGridEntity, id++, Traits::mydim);

        entityDofTuple.variableVector.assign(1, dofType);
        entityDofTuple.entityType = EntityType::vertex;
      }
      return entDofTupleVector;
    }

  private:
    GridElementEntityType elementGridEntity;
    IndexSetType const *const indexSet_;
  };
}  // namespace Ikarus::FiniteElements
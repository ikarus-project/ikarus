//
// Created by lex on 18/12/2021.
//

#pragma once
//#include <dune/functions/functionspacebases/basistags.hh>

#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/utils/concepts.h>

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

//    /** \brief Type of the Variables */
//    using VariableVectorType = typename FERequirementType::VariableType;
//
//    /** \brief Type of the DataVector */
//    using DataVectorType = typename FERequirementType::DataType;

    /** \brief Type of the Dofs / SolutionType
     * using NodalSolutionType = Displacement<ctype,worlddim>;*/

    /** \brief Type of the internal forces */
    using VectorType = Eigen::VectorXd;

    /** \brief Type of the stiffness matrix */
    using MatrixType = Eigen::MatrixXd;
  };


  template <typename LocalView>
  class FEVertexDisplacement {
  public:
    using RootBasis = typename LocalView::GlobalBasis;
    explicit FEVertexDisplacement(const LocalView& p_localView) : localView{p_localView} {
      static_assert(Ikarus::Concepts::PowerBasis<RootBasis>,
                    "You didn't pass a localview of a power basis to this method");
      static_assert(RootBasis::PreBasis::Node::CHILDREN == worlddim,
                    "The power basis children number does not coincide with the world space where the grid entity is "
                    "embedded into!");
      static_assert(
          Ikarus::Concepts::FlatIndexBasis<RootBasis>,
          "The Index merging strategy of the basis you passed has to be FlatLexicographic or FlatInterLeaved");
    }

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using DofTupleVectorType    = typename IFiniteElement::DofPairVectorType;
    using GridElementEntityType = typename LocalView::Element;
    using Traits                = FETraits<GridElementEntityType>;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = Traits::worlddim;

    [[nodiscard]] constexpr int dofSizeImpl() const { return localView.size(); }

    [[nodiscard]] DofTupleVectorType getEntityVariableTuple() const {
      const auto& fe = localView.tree().child(0).finiteElement();
      using namespace Ikarus::Variable;
      VariableTags dofType;
      if constexpr (worlddim == 3)
        dofType = VariableTags::displacement3d;
      else if constexpr (worlddim == 2)
        dofType = VariableTags::displacement2d;
      else if constexpr (worlddim == 1)
        dofType = VariableTags::displacement1d;
      DofTupleVectorType dofTupleVectorType;
      for (size_t i = 0; i < fe.size(); ++i) {
        dofTupleVectorType.emplace_back();
        auto& currentVariable = dofTupleVectorType.back();
        currentVariable.variableVector.push_back(dofType);
        // The zeros index of the multiIndex returns the flat index of the displacement degrees of freedom
        for (int j = 0; j < worlddim; ++j) {
          currentVariable.indices.push_back(localView.index((localView.tree().child(j).localIndex(i)))[0]);
        }
      }
      return dofTupleVectorType;
    }

    const GridElementEntityType& getEntity()
    {
      return localView.element();
    }

  private:
    LocalView localView;
  };
}  // namespace Ikarus::FiniteElements
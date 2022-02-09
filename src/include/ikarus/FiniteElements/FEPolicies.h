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

    /** \brief Type of the internal forces */
    using VectorType = Eigen::VectorXd;

    /** \brief Type of the stiffness matrix */
    using MatrixType = Eigen::MatrixXd;
  };


  template <typename LocalView>
  class FEDisplacement {
  public:
    using RootBasis = typename LocalView::GlobalBasis;
    using GlobalIndex = typename  LocalView::MultiIndex;
    explicit FEDisplacement(const LocalView& p_localView) : localView{p_localView} {
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
    using GridElementEntityType = typename LocalView::Element;
    using Traits                = FETraits<GridElementEntityType>;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = Traits::worlddim;

    [[nodiscard]] constexpr int dofSizeImpl() const { return localView.size(); }

    [[nodiscard]] std::vector<GlobalIndex> globalIndices() const {
      const auto& fe = localView.tree().child(0).finiteElement();
      std::vector<GlobalIndex> globalIndices;
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < worlddim; ++j) {
          globalIndices.push_back(localView.index((localView.tree().child(j).localIndex(i))));
        }
      }
      return globalIndices;
    }

    const GridElementEntityType& getEntity()
    {
      return localView.element();
    }

  private:
    LocalView const& localView;
  };
}  // namespace Ikarus::FiniteElements
//
// Created by lex on 18/12/2021.
//

#pragma once

#include "interfaceFiniteElement.hh"

#include <ikarus/utils/concepts.hh>

namespace Ikarus{

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

    /** \brief Type of the internal forces */
    using ScalarType = ctype;

    /** \brief Type of the stiffness matrix */
    using MatrixType = Eigen::MatrixXd;
  };




}  // namespace Ikarus::FiniteElements
// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/utils/concepts.hh>

namespace Ikarus {

  template <typename GridElementEntityType, bool useRef = false>
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
    using VectorType = std::conditional_t<useRef, Eigen::Ref<Eigen::VectorXd>, Eigen::VectorXd>;

    /** \brief Type of the internal forces */
    using ScalarType = ctype;

    /** \brief Type of the stiffness matrix */
    using MatrixType = std::conditional_t<useRef, Eigen::Ref<Eigen::MatrixXd>, Eigen::MatrixXd>;
  };

}  // namespace Ikarus

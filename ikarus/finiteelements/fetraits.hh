// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file fetraits.hh
 * \brief FETraits template structure for finite element traits.
 */

#pragma once

#include <ikarus/utils/concepts.hh>

namespace Ikarus {

  /**
   * \brief Template structure defining traits for a given grid element entity type.
   *
   * \tparam GridElement Type of the grid element entity.
   * \tparam useRef Boolean indicating whether to use Eigen::Ref for VectorType and MatrixType.
   */
  template <typename GridElement, bool useRef = false>
  struct FETraits {
    /** \brief Type used for coordinates */
    using ctype = double;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = GridElement::Geometry::coorddimension;

    /** \brief Dimension of the geometry */
    static constexpr int mydim = GridElement::mydimension;

    /** \brief Dimension of the grid */
    static constexpr int dimension = GridElement::dimension;

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

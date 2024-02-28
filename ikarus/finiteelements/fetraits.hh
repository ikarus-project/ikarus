// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file fetraits.hh
 * \brief FETraits template structure for finite element traits.
 */

#pragma once

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus {

/**
 * \brief Traits for handling finite elements.
 *
 * \tparam BH The basis handler type for the finite element.
 * \tparam useFlat A boolean indicating if the type of the underlying basis is of the flat or the untouched version.
 * \tparam FER The requirements for the finite element.
 * \tparam useRef Boolean indicating whether to use Eigen::Ref for VectorType and MatrixType.
 */
template <typename BH, bool useFlat, typename FER = FERequirements<>, bool useRef = false>
struct FETraits
{
  /** \brief Type of the basis handler of the finite element */
  using BasisHandler = BH;

  /** \brief Type of the requirements for the finite element */
  using FERequirementType = FER;

  /** \brief Type of the flat basis */
  using FlatBasis = typename BasisHandler::FlatBasis;

  /** \brief Type of the untouched basis */
  using UntouchedBasis = typename BasisHandler::UntouchedBasis;

  /** \brief Type of the basis version*/
  using Basis = std::conditional_t<useFlat, FlatBasis, UntouchedBasis>;

  /** \brief Type of the local view */
  using LocalView = typename Basis::LocalView;

  /** \brief Type of the grid view */
  using GridView = typename Basis::GridView;

  /** \brief Type of the grid element */
  using Element = typename LocalView::Element;

  /** \brief Type of the element geometry */
  using Geometry = typename Element::Geometry;

  /** \brief Type of the global index */
  using GlobalIndex = typename LocalView::MultiIndex;

  /** \brief Type used for coordinates */
  using ctype = double;

  /** \brief Dimension of the world space */
  static constexpr int worlddim = Geometry::coorddimension;

  /** \brief Dimension of the geometry */
  static constexpr int mydim = Element::mydimension;

  /** \brief Dimension of the grid */
  static constexpr int dimension = Element::dimension;

  /** \brief Type of the  coordinate */
  using GlobalCoordinates = Eigen::Matrix<ctype, worlddim, 1>;

  /** \brief Type of the ParameterSpace coordinate */
  using ParameterSpaceType = Eigen::Matrix<ctype, mydim, 1>;

  /** \brief Type of the internal forces */
  template <typename ScalarType = ctype>
  using VectorType = std::conditional_t<useRef, Eigen::Ref<Eigen::VectorX<ScalarType>>, Eigen::VectorX<ScalarType>&>;

  /** \brief Type of the stiffness matrix */
  template <typename ScalarType = ctype>
  using MatrixType = std::conditional_t<useRef, Eigen::Ref<Eigen::MatrixX<ScalarType>>, Eigen::MatrixX<ScalarType>&>;
};

} // namespace Ikarus

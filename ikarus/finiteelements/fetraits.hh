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
 * \tparam B The basis type for the finite element.
 * \tparam FER The requirements for the finite element.
 * \tparam useRef Boolean indicating whether to use Eigen::Ref for VectorType and MatrixType.
 */
template <typename B, typename FER = FERequirements<>, bool useRef = false>
struct FETraits
{
  /** \brief Type of the basis of the finite element */
  using Basis = B;

  /** \brief Type of the requirements for the finite element */
  using FERequirementType = FER;

  /** \brief Type of the flat basis */
  using FlatBasis = typename Basis::FlatBasis;

  /** \brief Type of the local view */
  using LocalView = typename FlatBasis::LocalView;

  /** \brief Type of the grid view */
  using GridView = typename FlatBasis::GridView;

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

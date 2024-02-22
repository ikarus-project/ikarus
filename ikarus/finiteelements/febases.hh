// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file febases.hh
 * \brief Contains the FEBase class, which works with different types of bases.
 */

#pragma once

#include <ikarus/finiteelements/fetraits.hh>

namespace Ikarus {

/**
 * \brief FEBase class for working with different types of bases.
 *
 * \details A single class which can be used by the finite elements to get different information
 * from the local view of the element.
 *
 * \tparam B The type of the basis.
 */
template <typename B>
class FEBase
{
public:
  using Traits      = FETraits<B>;                  ///< Type of the traits.
  using LocalView   = typename Traits::LocalView;   ///< Type of the local view.
  using GlobalIndex = typename Traits::GlobalIndex; ///< Type of the global index.
  using GridElement = typename Traits::Element;     ///< Type of the grid element.

  /**
   * \brief Constructor for the FEBase class.
   *
   * \param basis The basis.
   * \param element The local element.
   */
  explicit FEBase(const B& basis, const typename LocalView::Element& element)
      : localView_{basis.flat().localView()} {
    localView_.bind(element);
  }

  /**
   * \brief Get the size of the local view.
   * \return The size of the local view.
   */
  [[nodiscard]] constexpr size_t size() const { return localView_.size(); }

  /**
   * \brief Get the grid element associated with the local view.
   * \return The grid element.
   */
  const GridElement& gridElement() const { return localView_.element(); }

  /**
   * \brief Get the const reference to the local view.
   * \return The const reference to the local view.
   */
  const LocalView& localView() const { return localView_; }

  /**
   * \brief Get the reference to the local view.
   * \return The reference to the local view.
   */
  LocalView& localView() { return localView_; }

  /**
   * \brief Get the number of children of the basis.
   *
   * \return The number of children of the basis.
   */
  [[nodiscard]] int numberOfChildren() const { return localView_.tree().degree(); }

private:
  LocalView localView_;
};
} // namespace Ikarus

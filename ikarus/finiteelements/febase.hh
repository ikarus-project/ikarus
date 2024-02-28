// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file febase.hh
 * \brief Contains the FEBase class, which is used as a base class for all finite elements.
 * It provides information about the local view of a finite element.
 */

#pragma once

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/fetraits.hh>

namespace Ikarus {

/**
 * \brief FEBase class is a base class for all finite elements.
 *
 * \details A single class which can be used by the finite elements to get different information
 * from the local view of the element. It works for both the flat and untouched version of the basis.
 *
 * \tparam BH The type of the basis handler.
 * \tparam useFlat A boolean indicating if the type of the underlying basis is of the flat or the untouched version.
 * \tparam FER The requirements for the finite element.
 * \tparam useEigenRef A boolean flag indicating whether to use Eigen references.
 */
template <typename BH, bool useFlat = true, typename FER = FERequirements<>, bool useEigenRef = false>
class FEBase
{
public:
  using Traits      = FETraits<BH, useFlat, FER, useEigenRef>; ///< Type of the traits.
  using LocalView   = typename Traits::LocalView;              ///< Type of the local view.
  using GlobalIndex = typename Traits::GlobalIndex;            ///< Type of the global index.
  using GridElement = typename Traits::Element;                ///< Type of the grid element.

  /**
   * \brief Constructor for the FEBase class.
   *
   * \param basis The basis.
   * \param element The local element.
   */
  explicit FEBase(const BH& basisHandler, const typename LocalView::Element& element)
      : localView_{[&]() -> decltype(auto) {
          if constexpr (useFlat)
            return basisHandler.flat().localView();
          else
            return basisHandler.untouched().localView();
        }()} {
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

private:
  LocalView localView_;
};
} // namespace Ikarus

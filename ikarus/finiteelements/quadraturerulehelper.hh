// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file quadraturerulehelper.hh
 * \brief Helper functions to set quadrature rules for integration in finite elements.
 */

#pragma once

#include <dune/geometry/quadraturerules.hh>

namespace Ikarus {

/**
 * \brief A helper function to provide a default quadrature rule for integration.
 *
 * \details A helper function to provide a default (QuadratureType::GaussLegendre) quadrature rule for integration. It
 * also provides a default, if trimming operations are performed while working with NURBS basis functions.
 * \tparam dim Dimension of the grid element (mydimension)
 * \tparam Element Type of the grid element.
 * \param element The grid element.
 * \param order The order of functions to be integrated.
 * \return A default quadrature rule.
 */
template <int dim, typename Element>
auto defaultQuadratureRule(const Element& element, int order) {
  const auto defaultRule = Dune::QuadratureRules<double, dim>::rule(element.type(), order);
  if constexpr (requires { element.impl().getQuadratureRule(order); })
    if (element.impl().isTrimmed())
      return element.impl().getQuadratureRule(order);
    else
      return defaultRule;
  else
    return defaultRule;
}

/**
 * \brief A helper function to convert the number of Gauss integration points to the polynomial order of the function
 * being integrated.
 * \details Dune::QuadratureRules usually uses order of a function as an argument instead of number of
 * integration points (nGP). This function takes in nGP as an argument and returns the polynomial order (2 * nGP - 1)
 * upto which a function could be exactly integrated.
 * \param nGP Number of Gauss points to be used.
 * \return Polynomial order upto which a function can be exactly integrated.
 */
int numberOfGaussPointsToOrder(int nGP) { return 2 * nGP - 1; }

} // namespace Ikarus

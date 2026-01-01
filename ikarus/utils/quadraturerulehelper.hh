// SPDX-FileCopyrightText: 2021-2026 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file quadraturerulehelper.hh
 * \brief Free function for creating a tensor product quadrature rule and other Dune::Quadrature rule utils
 */
#pragma once

#include <dune/geometry/quadraturerules.hh>
namespace Ikarus {

/**
 * \brief Computes the tensor product quadrature rule using a base and one-dimensional quadrature rule.
 *
 * \details This function constructs a tensor product quadrature rule by combining a base quadrature rule
 * and a one-dimensional quadrature rule. It iterates over the points of the base quadrature and
 * combines them with each point of the one-dimensional quadrature.
 *
 * \tparam BaseQuadrature Type of the base quadrature.
 * \tparam Quadrature Type of the one-dimensional quadrature.
 *
 * \param baseQuad The base quadrature rule.
 * \param onedQuad The one-dimensional quadrature rule.
 *
 * \return Dune::QuadratureRule<double, BaseQuadrature::d + 1> - Tensor product quadrature rule.
 */
template <class BaseQuadrature, class Quadrature>
auto tensorProductQuadrature(const BaseQuadrature& baseQuad, const Quadrature& onedQuad) {
  constexpr int baseQuadDim       = BaseQuadrature::d;
  auto rule                       = Dune::QuadratureRule<double, baseQuadDim + 1>();
  const unsigned int baseQuadSize = baseQuad.size();
  for (unsigned int bqi = 0; bqi < baseQuadSize; ++bqi) {
    const typename Dune::QuadraturePoint<double, baseQuadDim>::Vector& basePoint = baseQuad[bqi].position();
    const double& baseWeight                                                     = baseQuad[bqi].weight();

    typename Dune::QuadraturePoint<double, baseQuadDim + 1>::Vector point;
    for (unsigned int i = 0; i < baseQuadDim; ++i)
      point[i] = basePoint[i];

    const unsigned int onedQuadSize = onedQuad.size();
    for (unsigned int oqi = 0; oqi < onedQuadSize; ++oqi) {
      point[baseQuadDim] = onedQuad[oqi].position()[0];
      rule.emplace_back(Dune::QuadraturePoint(point, baseWeight * onedQuad[oqi].weight()));
    }
  }
  return rule;
}

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

// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <dune/geometry/quadraturerules.hh>
namespace Ikarus
{
  template <class BaseQuadrature, class Quadrature>
  auto tensorProductQuadrature(const BaseQuadrature & baseQuad, const Quadrature & onedQuad)
  {
    constexpr int baseQuadDim= BaseQuadrature::d;
    auto rule = Dune::QuadratureRule<double,baseQuadDim+1>() ;
    const unsigned int baseQuadSize = baseQuad.size();
    for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
    {
      const typename Dune::QuadraturePoint<double, baseQuadDim>::Vector &   basePoint = baseQuad[bqi].position( );
      const double &baseWeight = baseQuad[bqi].weight( );

      typename Dune::QuadraturePoint<double, baseQuadDim+1>::Vector point;
      for( unsigned int i = 0; i < baseQuadDim; ++i )
        point[ i ] = basePoint[ i ];

      const unsigned int onedQuadSize = onedQuad.size();
      for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
      {
        point[ baseQuadDim ] = onedQuad[oqi].position()[ 0 ];
        rule.emplace_back( Dune::QuadraturePoint(point, baseWeight * onedQuad[oqi].weight()) );
      }
    }
    return rule;
  }


  auto transformFromZeroOneToMinusOneToOne(const Dune::QuadraturePoint<double,1>& gp)
  {
    typename Dune::QuadraturePoint<double,1>::Vector mappedCoordinate = gp.position();
    mappedCoordinate[0]= 2*mappedCoordinate[0]-1.0;
    double mappedWeight = gp.weight() * 2.0;
    return Dune::QuadraturePoint<double,1>(mappedCoordinate, mappedWeight);
  }
}

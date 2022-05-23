//
// Created by alex on 3/17/22.
//

#include "meta.hh"
namespace Ikarus::DerivativeDirections {
  SpatialPartial spatial(size_t i) { return {i}; }

  SingleCoeff<0> coeff(size_t i) {
    using namespace Dune::Indices;
    SingleCoeff<0> coeffs;
    std::get<1>(coeffs.index[_0]._data) = i;
    return coeffs;
  }

  TwoCoeff<0, 0> coeff(size_t i, size_t j) {
    using namespace Dune::Indices;
    TwoCoeff<0, 0> coeffs;
    std::get<1>(coeffs.index[_0]._data) = i;
    std::get<1>(coeffs.index[_1]._data) = j;
    return coeffs;
  }
}  // namespace Ikarus::DerivativeDirections

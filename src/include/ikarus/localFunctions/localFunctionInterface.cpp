//
// Created by alex on 3/17/22.
//

#include "localFunctionInterface.hh"
namespace Ikarus::DerivativeDirections {
  SpatialPartial spatial(size_t i) { return {i}; }
SingleCoeff coeff(size_t i) { return {i}; }
TwoCoeff coeff(size_t i, size_t j) { return {i,j}; }
}  // namespace Ikarus::DerivativeDirections

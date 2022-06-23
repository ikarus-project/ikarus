/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */



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

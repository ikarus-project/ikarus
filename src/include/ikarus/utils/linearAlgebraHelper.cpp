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

#include "linearAlgebraHelper.hh"

namespace Ikarus {
  Ikarus::DerivativeDirections::DerivativeNoOp transpose(const Ikarus::DerivativeDirections::DerivativeNoOp&) {
    return Ikarus::DerivativeDirections::DerivativeNoOp();
  }

  Ikarus::DerivativeDirections::DerivativeNoOp operator+(Ikarus::DerivativeDirections::DerivativeNoOp,
                                                         Ikarus::DerivativeDirections::DerivativeNoOp) {
    return Ikarus::DerivativeDirections::DerivativeNoOp();
  }

  Ikarus::DerivativeDirections::DerivativeNoOp operator-(Ikarus::DerivativeDirections::DerivativeNoOp,
                                                         Ikarus::DerivativeDirections::DerivativeNoOp) {
    return Ikarus::DerivativeDirections::DerivativeNoOp();
  }

  Ikarus::DerivativeDirections::DerivativeNoOp eval(const Ikarus::DerivativeDirections::DerivativeNoOp&) {
    return Ikarus::DerivativeDirections::DerivativeNoOp();
  }

  Ikarus::DerivativeDirections::DerivativeNoOp operator-(Ikarus::DerivativeDirections::DerivativeNoOp) {
    return Ikarus::DerivativeDirections::DerivativeNoOp();
  }
}  // namespace Ikarus
// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

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

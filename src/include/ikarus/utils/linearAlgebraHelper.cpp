//
// Created by Alex on 28.04.2022.
//

#include "linearAlgebraHelper.hh"

namespace Ikarus {
Ikarus::DerivativeDirections::DerivativeNoOp transpose(const Ikarus::DerivativeDirections::DerivativeNoOp&) {
  return Ikarus::DerivativeDirections::DerivativeNoOp();
}

Ikarus::DerivativeDirections::DerivativeNoOp operator +(Ikarus::DerivativeDirections::DerivativeNoOp,Ikarus::DerivativeDirections::DerivativeNoOp)
{
  return Ikarus::DerivativeDirections::DerivativeNoOp();
}

Ikarus::DerivativeDirections::DerivativeNoOp operator -(Ikarus::DerivativeDirections::DerivativeNoOp,Ikarus::DerivativeDirections::DerivativeNoOp)
{
  return Ikarus::DerivativeDirections::DerivativeNoOp();
}

Ikarus::DerivativeDirections::DerivativeNoOp eval(const Ikarus::DerivativeDirections::DerivativeNoOp&)
{
  return Ikarus::DerivativeDirections::DerivativeNoOp();
}
}
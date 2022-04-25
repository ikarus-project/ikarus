//
// Created by Alex on 21.04.2022.
//

#pragma once

#include <ikarus/localFunctions/localFunctionInterface.hh>
namespace Ikarus
{

template <typename Op, typename E1,typename E2>
class BinaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op> {

 protected:
  E1 const& _u;
  E2 const& _v;
 public:
  BinaryLocalFunctionExpression(Ikarus::LocalFunctionInterface<E1> const& u, Ikarus::LocalFunctionInterface<E2> const& v) : _u(static_cast<E1 const&>(u)), _v(static_cast<E2 const&>(v)) {  }

  static constexpr bool isLeaf = false;
  template<int i = 0>
  const auto& basis()const
  {
    if constexpr(i==0)
      return _u.basis();
    else
      return _v.basis();
  }

 protected:
  template<int i, typename LocalFunctionEvaluationArgs_>
  auto evalDerivative(const LocalFunctionEvaluationArgs_& localFunctionArgs) const
  {
    if constexpr(i==0)
      return evaluateDerivativeImpl(_u,localFunctionArgs);
    else
      return evaluateDerivativeImpl(_v,localFunctionArgs);
  }


};

}
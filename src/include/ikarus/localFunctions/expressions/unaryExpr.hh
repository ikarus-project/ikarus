//
// Created by Alex on 21.04.2022.
//

#pragma once

#include <ikarus/localFunctions/localFunctionInterface.hh>
namespace Ikarus
{

template <typename Op, typename E1>
class UnaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op> {

 protected:
  E1 const& m;
 public:
   UnaryLocalFunctionExpression(Ikarus::LocalFunctionInterface<E1> const& u) : m(static_cast<E1 const&>(u)) {  }

  static constexpr bool isLeaf = false;

};

}
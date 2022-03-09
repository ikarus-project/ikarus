//
// Created by lex on 08/03/2022.
//

#pragma once

#include <dune/functions/analyticfunctions/polynomial.hh>

#include "Eigen/Dense"
namespace Ikarus {
  std::tuple<Dune::Functions::Polynomial<double>,double> polyfit(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>  &y, const int deg);
}  // namespace Ikarus
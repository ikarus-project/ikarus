//
// Created by lex on 08/03/2022.
//

#pragma once
#include "polyfit.hh"

#include <Eigen/Core>

namespace Ikarus {
  /*
   * This function is inspired from
   * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/identify_linear_piece.m
   */
  std::tuple<Dune::Functions::Polynomial<double>, decltype(Eigen::seq(0,0))>
  findLineSegment(const Eigen::VectorXd& x, const Eigen::VectorXd& y, int segmentSize);
}  // namespace Ikarus
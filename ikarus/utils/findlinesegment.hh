// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "polyfit.hh"

#include <Eigen/Core>

namespace Ikarus {
  /*
   * This function is inspired from
   * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/identify_linear_piece.m
   */
  std::tuple<Dune::Functions::Polynomial<double>, decltype(Eigen::seq(0, 0))> findLineSegment(const Eigen::VectorXd& x,
                                                                                              const Eigen::VectorXd& y,
                                                                                              int segmentSize);
}  // namespace Ikarus

// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/functions/analyticfunctions/polynomial.hh>

#include <Eigen/Dense>
namespace Ikarus {
  std::tuple<Dune::Functions::Polynomial<double>, double> polyfit(const Eigen::Ref<const Eigen::VectorXd>& x,
                                                                  const Eigen::Ref<const Eigen::VectorXd>& y,
                                                                  const int deg);
}  // namespace Ikarus

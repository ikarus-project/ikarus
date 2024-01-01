// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>

namespace Ikarus {

  struct ControlInformation {
    bool success{false};
    std::vector<Ikarus::NonLinearSolverInformation> solverInfos{};
    int totalIterations{0};
  };
}  // namespace Ikarus

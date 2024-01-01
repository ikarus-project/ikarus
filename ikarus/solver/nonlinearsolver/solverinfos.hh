// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace Ikarus {
  struct NonLinearSolverInformation {
    explicit operator bool() const { return success; }
    bool success{false};
    double residualNorm{std::numeric_limits<double>::infinity()};
    double correctionNorm{std::numeric_limits<double>::infinity()};
    int iterations{-1};
  };
}  // namespace Ikarus

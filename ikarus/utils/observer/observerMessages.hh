// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

enum class ControlMessages { BEGIN, CONTROL_STARTED, CONTROL_ENDED, STEP_STARTED, STEP_ENDED, SOLUTION_CHANGED, END };

enum class NonLinearSolverMessages {
  BEGIN,
  INIT,
  ITERATION_STARTED,
  ITERATION_ENDED,
  RESIDUALNORM_UPDATED,
  CORRECTIONNORM_UPDATED,
  SCALARSUBSIDIARY_UPDATED,
  SOLUTION_CHANGED,
  FINISHED_SUCESSFULLY,
  END
};

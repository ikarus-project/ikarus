//
// Created by alex on 12/20/21.
//

#pragma once

enum class ControlMessages {
  BEGIN,
  CONTROL_STARTED,
  CONTROL_ENDED,
  LOADSTEP_STARTED,
  LOADSTEP_ENDED,
   SOLUTION_CHANGED,
  END
};

enum class NonLinearSolverMessages {
  BEGIN,
  ITERATION_STARTED,
  ITERATION_ENDED,
  RESIDUALNORM_UPDATED,
  SOLUTION_CHANGED,
  END
};
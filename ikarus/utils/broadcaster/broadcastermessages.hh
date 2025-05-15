// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file broadcastermessages.hh
 * \brief Enums for observer messages
 * \ingroup observer
 */

#pragma once
#include <ikarus/utils/makeenum.hh>

namespace Ikarus {

/**
 * \brief Enum class defining control-routine-related messages.
 */
MAKE_ENUM(ControlMessages, CONTROL_STARTED, CONTROL_ENDED, STEP_STARTED, STEP_ENDED, SOLUTION_CHANGED)

/**
 * \brief Enum class defining non-linear solver-related messages.
 */
MAKE_ENUM(NonLinearSolverMessages, INIT, ITERATION_STARTED, ITERATION_ENDED, RESIDUALNORM_UPDATED,
          CORRECTIONNORM_UPDATED, CORRECTION_UPDATED, SOLUTION_CHANGED, FINISHED_SUCESSFULLY);

} // namespace Ikarus

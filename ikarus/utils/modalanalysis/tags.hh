// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file tags.hh
 * \brief Tags for ModalAnalysis related things
 */

#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

/**
 * \brief A strongly typed enum class representing the type of result of a modal analysis
 */
MAKE_ENUM(ModalAnalysisResultType, squaredAngularFrequency, angularFrequency, naturalFrequency);

} // namespace Ikarus::Dynamics
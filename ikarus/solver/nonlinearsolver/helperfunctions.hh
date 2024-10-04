// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file helperfunctions.hh
 * \brief Implementation of the helper functions used by the nonlinear solvers.
 */

#pragma once

#include <ikarus/utils/concepts.hh>
namespace Ikarus {
template <typename Assembler>
void updateStates(std::shared_ptr<Assembler>& assembler, const auto& correction) {
  if constexpr (Concepts::FlatAssembler<Assembler>) {
    auto fes = assembler->finiteElements();
    auto req = assembler->requirement();
    for (auto& fe : fes)
      updateState(fe, req, correction);
  } else {
    /* Dummy branch used is std::is_same_v<Assembler, Impl::NoAssembler> */
  }
}
} // namespace Ikarus
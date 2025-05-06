// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlroutinestate.hh
 * \brief State for all controlroutines
 */

#pragma once

#include <ikarus/utils/traits.hh>

namespace Ikarus {

/**
 * \brief State for control routines
 *
 * \tparam D the type of the domain (in most cases FERequirement)
 */
template <typename D>
struct ControlRoutineState
{
  using Domain = D;

  const Domain& domain;
  int loadStep{};
  double stepSize{};
};

namespace Impl {

  template <typename F>
  struct ControlRoutineStateFactory
  {
  private:
    using SignatureTraits = typename F::Traits;
    using Domain          = typename SignatureTraits::Domain;

  public:
    using type = ControlRoutineState<Domain>;
  };
} // namespace Impl

/**
 * \brief Helper to deduce the correct types for ControlRoutineState
 *
 * \tparam F Type of the differentiable function to solve.
 */
template <typename F>
using ControlRoutineStateType = Impl::ControlRoutineStateFactory<F>::type;

} // namespace Ikarus

// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverstate.hh
 * \brief State for all nonlinear solvers
 */

#pragma once

#include <ikarus/utils/traits.hh>

namespace Ikarus {

/**
 * \brief State for nonlinear solvers
 *
 * \tparam D the type of the domain (in most cases FERequirement)
 * \tparam CT the type of the correction vector
 */
template <typename D, typename CT>
struct NonlinearSolverState
{
  using Domain         = D;
  using CorrectionType = CT;

  const Domain& domain;
  const CorrectionType& correction;

  double rNorm{};
  double dNorm{};
  int iteration{};
};

namespace Impl {

  template <typename SignatureTraits, int n>
  struct CorrectionType;

  template <typename SignatureTraits>
  struct CorrectionType<SignatureTraits, 1>
  {
    using type = typename SignatureTraits::template Range<0>;
  };

  template <typename SignatureTraits>
  struct CorrectionType<SignatureTraits, 2>
  {
    using type = typename SignatureTraits::template Range<1>;
  };

  template <typename F>
  struct NonlinearSolverStateFactory
  {
  private:
    using SignatureTraits = typename F::Traits;
    using Domain          = typename SignatureTraits::Domain;

    constexpr static int nDerivatives = F::nDerivatives;

  public:
    using type = NonlinearSolverState<Domain, typename CorrectionType<SignatureTraits, nDerivatives>::type>;
  };
} // namespace Impl

/**
 * \brief Helper to deduce the correct types for NonlinearSolverState
 *
 * \tparam F Type of the differentiable function to solve.
 */
template <typename F>
using NonlinearSolverStateType = Impl::NonlinearSolverStateFactory<F>::type;

} // namespace Ikarus

// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverstate.hh
 * \brief State for all nonlinear solvers
 */

#pragma once

#include <ikarus/utils/traits.hh>

namespace Ikarus {

// Forward declarations
#ifndef DOXYGEN
template <typename TypeListOne, typename TypeListTwo>
class NonLinearOperator;
#endif

/**
 * \brief State for nonlinear solvers
 *
 * \tparam Correction the type of the correction vector
 * \tparam SolutionType the type of the solution vector
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

  template <typename F>
  struct NonlinearSolverStateFactory
  {
  private:
    using SignatureTraits = typename F::Traits;
    using Domain          = typename SignatureTraits::Domain;

    template <int n>
    struct CorrectionType;

    template <>
    struct CorrectionType<3>
    {
      using type = typename SignatureTraits::template Range<0>;
    };

    template <>
    struct CorrectionType<4>
    {
      using type = typename SignatureTraits::template Range<1>;
    };

    constexpr static int numRanges = SignatureTraits::numberOfRanges;
    // static_assert(numRanges == 3);
    // using CorrectionType           = std::conditional_t<numRanges == 3, typename SignatureTraits::template Range<0>,
    //                                                     typename SignatureTraits::template Range<1>>;

  public:
    using type = NonlinearSolverState<Domain, typename CorrectionType<numRanges>::type>;
  };
} // namespace Impl

/**
 * \brief Helper to deduce the correct types for NonlinearSolverState
 *
 * \tparam NLO The nonlinear operator
 */
template <typename F>
using NonlinearSolverStateType = Impl::NonlinearSolverStateFactory<F>::type;

} // namespace Ikarus

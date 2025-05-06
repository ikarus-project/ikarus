// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file lambertw.hh
 * \brief Implementation of the zeroth branch of the lambertW function
 */

#pragma once
#include <cmath>
#include <limits>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <ikarus/utils/concepts.hh>

namespace Ikarus::util {

namespace Impl {
  template <typename T>
  T log1p(T x) {
    return log(T(1.0) + x);
  }
} // namespace Impl

/**
 * \brief Implementation of the principal branch of the Lambert-W function (branch 0 in the domain
 * \f$ [- \dfrac{1}{e}, \inf) \f$ ), that is defined as the inverse function of \f$ f: x \mapsto xe^x \f$.
 *  It is defined for inputs of \f$ [- \dfrac{1}{e}, \inf) \f$
 *
 * \details The implementation uses Halley's iterative method and is inspired by
 * https://github.com/JuliaMath/LambertW.jl/blob/master/src/LambertW.jl (licensed under MIT "Expat" License).
 *
 * \tparam ST the ScalarType (defaults to double).
 * \param z the input value.
 * \param maxIterations optionally define a maximmum of iteration.
 * \param eps optionally define the epsilon of the iteration process (detuls to machine epsilon).
 * \return ST the output value.
 */
template <typename ST = double>
ST lambertW0(ST z, int maxIterations = 20, ST eps = std::numeric_limits<ST>::epsilon()) {
  if constexpr (not Concepts::AutodiffScalar<ST>) {
    if (std::isnan(z))
      return std::numeric_limits<ST>::quiet_NaN();
    if (std::isinf(z))
      return z;
  }

  const ST branchPoint = -1.0 / std::exp(1.0);

  // If z equals -1/e then W(z) = -1.
  if (Dune::FloatCmp::eq(z, branchPoint))
    return -1.0;

  // For branch 0 the domain is z >= -1/e.
  if (z < branchPoint)
    DUNE_THROW(Dune::InvalidStateException, "lambertW0: z must be >= -1/e for branch 0");

  // Choose an initial guess x0. See https://en.wikipedia.org/wiki/Lambert_W_function.
  ST x0;
  if (Dune::FloatCmp::gt(z, ST(1.0))) {
    ST lx  = log(z);
    ST llx = log(lx);
    x0     = lx - llx - 0.5 * Impl::log1p<ST>(-llx / lx);
  } else {
    x0 = 0.567 * z;
  }

  // Begin Halley's iterative method. See https://en.wikipedia.org/wiki/Halley%27s_method.
  ST x        = x0;
  ST lastDiff = 0.0;

  for (int iter = 0; iter < maxIterations; ++iter) {
    const ST ex          = exp(x);
    const ST f           = x * ex - z;
    const ST fPrime      = ex * (x + 1.0);
    const ST fPrimePrime = ex * (x + 2.0);
    const ST denom       = fPrime - ((1.0 / 2.0) * (fPrimePrime / fPrime) * f);

    const ST newX = x - f / denom;
    const ST diff = abs(newX - x);

    // Check for convergence:
    if (Dune::FloatCmp::le<ST>(diff, 3 * eps * abs(x)) || Dune::FloatCmp::eq<ST>(diff, lastDiff))
      return newX;

    lastDiff = diff;
    x        = newX;
  }

  DUNE_THROW(Dune::MathError, "lambertW0: failed to converge within the maximum number of iterations");
}
} // namespace Ikarus::util

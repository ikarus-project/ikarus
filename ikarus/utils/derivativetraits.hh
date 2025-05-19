// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file derivativetraits.hh
 * \brief Contains derivative traits for common vbalue and derivative relations used for
 * makeDifferentiableFunctionFromCallables
 */

#pragma once

#include <dune/functions/common/differentiablefunctionfromcallables.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ikarus/utils/traits.hh>

namespace Ikarus {

namespace Impl {
  template <typename... Args>
  struct Functions;
} // namespace Impl

#ifndef DOXYGEN

enum class FEParameter;
enum class FESolutions;
template <FESolutions sol, FEParameter para, typename SV, typename PM>
class FERequirements;

/**
 * \brief Represents the traits obtained from a set of callables
 * \ingroup utils
 * \tparam TypeListOne The type list for the first set of functions.
 * \tparam TypeListTwo The type list for the second set of functions.
 */
template <typename TypeListOne, typename TypeListTwo>
class DerivativeTraitsFromCallables
{
public:
  DerivativeTraitsFromCallables(const TypeListOne&, const TypeListTwo&) {
    static_assert(!sizeof(TypeListOne),
                  "This type should not be instantiated. check that your arguments satisfies the template below");
  }
};

#endif

template <typename... DerivativeArgs, typename Arg>
struct DerivativeTraitsFromCallables<Impl::Functions<DerivativeArgs...>, Arg>
{
  /**
   * \brief Constructor for DerivativeTraitsFromCallables.
   *
   * \param derivativesFunctions The Functions object for derivative arguments.
   * \param parameterI The Parameter object for parameter arguments.
   */
  DerivativeTraitsFromCallables(const Impl::Functions<DerivativeArgs...>& derivativesFunctions, const Arg& parameterI) {
  }

  using Ranges    = std::tuple<std::invoke_result_t<DerivativeArgs, const Arg&>..., Dune::Functions::InvalidRange>;
  using RawRanges = std::tuple<std::remove_cvref_t<std::invoke_result_t<DerivativeArgs, const Arg&>>...,
                               Dune::Functions::InvalidRange>;

  using Domain     = std::remove_cvref_t<Arg>;
  using Signatures = std::tuple<std::invoke_result_t<DerivativeArgs, const Arg&>(const Arg&)...,
                                Dune::Functions::InvalidRange(const Arg&)>;
  using RawSignatures =
      std::tuple<std::remove_cvref_t<std::invoke_result_t<DerivativeArgs, const Arg&>>(std::remove_cvref_t<Arg>)...,
                 Dune::Functions::InvalidRange(std::remove_cvref_t<Arg>)>;

  template <int I>
  using Signature = std::tuple_element_t<I, Signatures>;

  template <int I>
  using Range = std::tuple_element_t<I, Ranges>;

  template <int I>
  using RawRange = std::tuple_element_t<I, RawRanges>;

  template <typename Signature>
  struct DerivativeTraits
  {
  private:
    static constexpr int indexOfSignature = traits::Index<Signature, RawSignatures>::value + 1;

  public:
    using Range = std::tuple_element_t<indexOfSignature, Ranges>;
  };
};

template <typename... DerivativeArgs, typename Arg>
DerivativeTraitsFromCallables(const Impl::Functions<DerivativeArgs...>& derivativesFunctions, Arg&& parameterI)
    -> DerivativeTraitsFromCallables<Impl::Functions<DerivativeArgs...>, Arg>;
} // namespace Ikarus

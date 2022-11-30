// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "localFunctionInterface.hh"
namespace Ikarus {

  /** Helper to evaluate the local basis ansatz function and gradient with an integration point index or coordinate
   * vector*/
  template <typename DomainTypeOrIntegrationPointIndex, typename Basis>
  auto evaluateFunctionAndDerivativeWithIPorCoord(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                  const Basis& basis) {
    if constexpr (std::is_same_v<DomainTypeOrIntegrationPointIndex, typename Basis::DomainType>) {
      typename Basis::JacobianType dN;
      basis.evaluateJacobian(localOrIpId, dN);
      typename Basis::AnsatzFunctionType N;
      basis.evaluateFunction(localOrIpId, N);
      return std::make_tuple(N, dN);
    } else if constexpr (std::numeric_limits<DomainTypeOrIntegrationPointIndex>::is_integer) {
      const typename Basis::JacobianType& dN      = basis.evaluateJacobian(localOrIpId);
      const typename Basis::AnsatzFunctionType& N = basis.evaluateFunction(localOrIpId);
      return std::make_tuple(std::ref(N), std::ref(dN));
    } else
      static_assert(
          std::is_same_v<DomainTypeOrIntegrationPointIndex,
                         typename Basis::DomainType> or std::is_same_v<DomainTypeOrIntegrationPointIndex, int>,
          "The argument you passed should be an id for the integration point or the point where the "
          "derivative should be evaluated");
  }

  /** Helper to evaluate the local basis ansatz function gradient with an integration point index or coordinate vector*/
  template <typename DomainTypeOrIntegrationPointIndex, typename Basis>
  auto evaluateDerivativeWithIPorCoord(const DomainTypeOrIntegrationPointIndex& localOrIpId, const Basis& basis) {
    if constexpr (std::is_same_v<DomainTypeOrIntegrationPointIndex, typename Basis::DomainType>) {
      typename Basis::JacobianType dN;
      basis.evaluateJacobian(localOrIpId, dN);
      return dN;
    } else if constexpr (std::numeric_limits<DomainTypeOrIntegrationPointIndex>::is_integer) {
      const typename Basis::JacobianType& dN = basis.evaluateJacobian(localOrIpId);
      return dN;
    } else
      static_assert(
          std::is_same_v<DomainTypeOrIntegrationPointIndex,
                         typename Basis::DomainType> or std::is_same_v<DomainTypeOrIntegrationPointIndex, int>,
          "The argument you passed should be an id for the integration point or the point where the "
          "derivative should be evaluated");
  }

  /** Helper to evaluate the local basis ansatz function with an integration point index or coordinate vector*/
  template <typename DomainTypeOrIntegrationPointIndex, typename Basis>
  auto evaluateFunctionWithIPorCoord(const DomainTypeOrIntegrationPointIndex& localOrIpId, const Basis& basis) {
    if constexpr (std::is_same_v<DomainTypeOrIntegrationPointIndex, typename Basis::DomainType>) {
      typename Basis::AnsatzFunctionType N;
      basis.evaluateFunction(localOrIpId, N);
      return N;
    } else if constexpr (std::numeric_limits<DomainTypeOrIntegrationPointIndex>::is_integer) {
      return basis.evaluateFunction(localOrIpId);
    } else
      static_assert(
          std::is_same_v<DomainTypeOrIntegrationPointIndex,
                         typename Basis::DomainType> or std::is_same_v<DomainTypeOrIntegrationPointIndex, int>,
          "The argument you passed should be an id for the integration point or the point where the "
          "derivative should be evaluated");
  }

  /** Helper to transform the derivatives if the transform arg is non-empty*/
  template <typename... TransformArgs>
  void maytransformDerivatives(const auto& dNraw, auto& dNTransformed,
                               const TransformWith<TransformArgs...>& transArgs) {
    if constexpr (sizeof...(TransformArgs) > 0) {
      dNTransformed = dNraw * std::get<0>(transArgs.args);
    } else
      dNTransformed = dNraw;
  }

}  // namespace Ikarus

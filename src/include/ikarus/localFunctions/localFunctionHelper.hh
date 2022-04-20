//
// Created by Alex on 20.04.2022.
//

#pragma once

namespace Ikarus
{
template <typename DomainTypeOrIntegrationPointIndex, typename Basis>
auto evaluateFunctionAndDerivativeWithIPorCoord(const DomainTypeOrIntegrationPointIndex& localOrIpId,const Basis& basis)  {
  if constexpr (std::is_same_v<DomainTypeOrIntegrationPointIndex, typename Basis::DomainType>) {
    typename Basis::JacobianType dN;
    basis.evaluateJacobian(localOrIpId, dN);
    typename Basis::AnsatzFunctionType N;
    basis.evaluateFunction(localOrIpId, N);
    return std::make_tuple(N, dN);
  } else if constexpr (std::numeric_limits<DomainTypeOrIntegrationPointIndex>::is_integer) {
    const typename Basis::JacobianType& dN = basis.evaluateJacobian(localOrIpId);
    const typename Basis::AnsatzFunctionType& N      = basis.evaluateFunction(localOrIpId);
    return std::make_tuple(std::ref(N), std::ref(dN));
  } else
    static_assert(std::is_same_v<DomainTypeOrIntegrationPointIndex,
                                 typename Basis::DomainType> or std::is_same_v<DomainTypeOrIntegrationPointIndex, int>,
                  "The argument you passed should be an id for the integration point or the point where the "
                  "derivative should be evaluated");
}

template <typename DomainTypeOrIntegrationPointIndex, typename Basis>
auto evaluateDerivativeWithIPorCoord(const DomainTypeOrIntegrationPointIndex& localOrIpId,const Basis& basis)  {
  if constexpr (std::is_same_v<DomainTypeOrIntegrationPointIndex, typename Basis::DomainType>) {
    typename Basis::JacobianType dN;
    basis.evaluateJacobian(localOrIpId, dN);
    return dN;
  } else if constexpr (std::numeric_limits<DomainTypeOrIntegrationPointIndex>::is_integer) {
    const typename Basis::JacobianType& dN = basis.evaluateJacobian(localOrIpId);
    return dN;
  } else
    static_assert(std::is_same_v<DomainTypeOrIntegrationPointIndex,
                                 typename Basis::DomainType> or std::is_same_v<DomainTypeOrIntegrationPointIndex, int>,
                  "The argument you passed should be an id for the integration point or the point where the "
                  "derivative should be evaluated");
}

template <typename DomainTypeOrIntegrationPointIndex, typename Basis>
auto evaluateFunctionWithIPorCoord(const DomainTypeOrIntegrationPointIndex& localOrIpId,const Basis& basis)  {
  if constexpr (std::is_same_v<DomainTypeOrIntegrationPointIndex, typename Basis::DomainType>) {
    typename Basis::AnsatzFunctionType N;
    basis.evaluateFunction(localOrIpId, N);
    return N;
  } else if constexpr (std::numeric_limits<DomainTypeOrIntegrationPointIndex>::is_integer) {
    return  basis.evaluateFunction(localOrIpId);
  } else
    static_assert(std::is_same_v<DomainTypeOrIntegrationPointIndex,
                                 typename Basis::DomainType> or std::is_same_v<DomainTypeOrIntegrationPointIndex, int>,
                  "The argument you passed should be an id for the integration point or the point where the "
                  "derivative should be evaluated");
}

template <typename... TransformArgs>
void maytransformDerivatives(const auto& dNraw, auto& dNTransformed,
                             const TransformWith<TransformArgs...>& transArgs)  {
  if constexpr (sizeof...(TransformArgs) > 0)
    dNTransformed = dNraw * std::get<0>(transArgs.args);
  else
    dNTransformed = dNraw;
}

}
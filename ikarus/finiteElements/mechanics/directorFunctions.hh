// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <Eigen/Core>
#include <dune/common/fvector.hh>
//#include <dune/localfefunctions/impl/standardLocalFunction.hh>
//#include <dune/localfefunctions/impl/projectionBasedLocalFunction.hh>
#include "NFE.hh"
#include "PBFE.hh"
#include "GFE.hh"
#include <dune/common/overloadset.hh>
#include <ranges>
#include <dune/localfefunctions/concepts.hh>
#include <variant>

namespace Ikarus {

template<class... Implementations>
class DirectorFunctionVariant
{
  using FirstImpl= std::tuple_element_t<0, std::tuple<Implementations...>>;

 public:

  template<class Implementation>
  explicit DirectorFunctionVariant(const Implementation& impl) :
      impl_(impl)
  {}

  DirectorFunctionVariant() = default;
  DirectorFunctionVariant(const DirectorFunctionVariant& other) = default;
  template<class Implementation> requires (!std::is_same_v<Implementation,DirectorFunctionVariant>)
  DirectorFunctionVariant& operator=(const Implementation& impl)
      {
        impl_=impl;
        return *this;
      };
  DirectorFunctionVariant(DirectorFunctionVariant&& other)  noexcept = default;
  DirectorFunctionVariant& operator=(const DirectorFunctionVariant& other) = default;
  DirectorFunctionVariant& operator=(DirectorFunctionVariant&& other)  noexcept = default;

  template <typename ...Args>
  auto evaluate(Args&&... args) const{

    return std::visit([&](const auto& impl) { return impl.evaluate(std::forward<Args>(args)...); }, impl_);
  }

  template <typename ...Args>
  auto evaluateDerivative(Args&&... args) const {

    return std::visit([&](const auto& impl) { return impl.evaluateDerivative(std::forward<Args>(args)...); }, impl_);
  }

//  template <typename... WrtArgs, typename Transform = Dune::DerivativeDirections::GridElement,
//            typename DomainTypeOrIntegrationPointIndex>
//  auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId, Dune::Wrt<WrtArgs...>&& args,
//                          const Dune::On<Transform>& transform = {}) const {
//    return std::visit([&](const auto& impl) { return impl.evaluateDerivative(localOrIpId,args,Dune::along(),transform); }, impl_);
//  }

//  template <typename IntegrationRule>
//  void bind(IntegrationRule&& p_rule, std::set<int>&& ints) {
//     std::visit([&](const auto& impl) { impl.basis.bind(std::forward<IntegrationRule>(p_rule), std::forward<std::set<int>>(ints)); }, impl_);
//  }

 private:
  std::variant<Implementations...> impl_;
};


template <typename DuneBasis, typename CoeffContainer, typename Geometry, std::size_t ID>
using DirectorFunctionVar = DirectorFunctionVariant< Dune::EmbeddedLocalFunction<DuneBasis,CoeffContainer,Geometry,ID>,
                                                        Dune::ProjectionBasedLocalFunction2<DuneBasis,CoeffContainer,Geometry,ID>,
                                                        Dune::GeodesicLocalFunction<DuneBasis,CoeffContainer,Geometry,ID>
                                                    >;

}

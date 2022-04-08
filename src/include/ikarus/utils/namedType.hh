//
// Created by lex on 15/03/2022.
//

#pragma once
#include <memory>

namespace Ikarus {
  /*
   * See https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
   */
  template <typename Parameter, typename... T>
  class NamedType {
    using CommonType = std::common_type_t<T...>;
    using ValueType  = std::conditional_t<sizeof...(T) == 1, CommonType, std::array<CommonType, sizeof...(T)>>;

  public:
    explicit NamedType(ValueType const& value) : value_(value) {}
    explicit NamedType(ValueType&& value) : value_(std::move(value)) {}
    ValueType& get() { return value_; }
    ValueType const& get() const { return value_; }

  private:
    ValueType value_;
  };
}  // namespace Ikarus
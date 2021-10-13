//
// Created by Alex on 27.05.2021.
//

#pragma once

namespace Ikarus::utils {
  /** \brief A traits which returns false if template is instiantiated, handy für static_assert */
  template <typename... Args>
  bool dependentFalse() {
    return false;
  }

  template <typename T>
  struct is_std_array : std::false_type {};

  template <typename T, std::size_t N>
  struct is_std_array<std::array<T, N>> : std::true_type {};
}  // namespace Ikarus::utils
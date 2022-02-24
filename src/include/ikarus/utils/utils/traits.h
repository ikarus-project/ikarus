//
// Created by Alex on 27.05.2021.
//

#pragma once

namespace Ikarus::utils {

  template <typename T>
  struct is_std_array : std::false_type {};

  template <typename T, std::size_t N>
  struct is_std_array<std::array<T, N>> : std::true_type {};
}  // namespace Ikarus::utils
// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file integersequence.hh
 * \brief integersequence helper
 */

#pragma once

#include <type_traits>
#include <utility>
#include <array>

#include <dune/common/indices.hh>
#include <dune/common/integersequence.hh>


// we put this function into the dune namespace, since in <dune/common/integersequence.hh> similar function reside
namespace Dune {
template <typename T, T... I>
consteval auto makeArrayFromIndexSequence(std::integer_sequence<T, I...> ={}) {
    return std::array<T, sizeof...(I)>{I...};
}
  namespace Impl{

template <auto firstSequence, auto secondSequence>
consteval auto set_intersectionImpl() {
    auto in1 = makeArrayFromIndexSequence(firstSequence);
    auto in2 = makeArrayFromIndexSequence(secondSequence);
    std::array<int, std::min(in1.size(), in2.size())> out{};
    auto first1 = in1.begin();
    auto last1 = in1.end();
    auto first2 = in2.begin();
    auto last2 = in2.end();
    auto d_first = out.begin();
    int counterEqual = 0;
    // set_intersection implementation
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2)
            ++first1;
        else {
            if (!(*first2 < *first1)) {
                *d_first++ = *first1++;  // *first1 and *first2 are equivalent.
                ++counterEqual;
            }
            ++first2;
        }
    }

    return std::make_pair(out, counterEqual);
}

}

template <class T, T... I, T... J>
constexpr auto set_intersection(std::integer_sequence<T, I...> seq1,
                                std::integer_sequence<T, J...> seq2) {
    constexpr auto pair = Impl::set_intersectionImpl<seq1, seq2>();

    auto iSeq = std::make_index_sequence<pair.second>{};
    return Dune::unpackIntegerSequence(
        [&](auto... i) { return std::integer_sequence<T, pair.first[i]...>{}; },
        iSeq);
}
} // namespace Dune

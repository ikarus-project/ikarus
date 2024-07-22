// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <ikarus/utils/integersequence.hh>


template<typename T, T... I,T... J,T... K>
constexpr void testIntegerSequence(std::integer_sequence<T,I...> a, std::integer_sequence<T,J...> b,std::integer_sequence<T,K...> expectedResult) {

    auto out = Dune::set_intersection(a, b);
    static_assert(expectedResult.size()==out.size());
    static_assert(Dune::equal(out,expectedResult));
}

  template<size_t... I>
  using seq = std::integer_sequence<size_t, I...>;

int main()
{

  testIntegerSequence(seq<1,2,3>{},seq<3,4,5>{},seq<3>{});
  testIntegerSequence(seq<1,3>{},seq<3,4,5>{},seq<3>{});
  testIntegerSequence(seq<1,3,5,17>{},seq<3,4,5,7,9,17>{},seq<3,5,17>{});

}

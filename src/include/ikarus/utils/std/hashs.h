//
// Created by Alex on 26.06.2021.
//

#pragma once
#include <concepts>
/*
 *  This is an implemation of a has functions for std::pair
 *  It uses the XOR operator to combine the hashes
 *  This is only good for distinct types T1,T2 since it avoid collisions
 *  https://www.techiedelight.com/use-std-pair-key-std-unordered_map-cpp/
 *  https://stackoverflow.com/questions/39370214/how-to-properly-hash-a-pair-of-pointers
 */
struct pair_hash_naive
{
  template <class T1, class T2> requires (!std::same_as<T1,T2>)
  std::size_t operator() (const std::pair<T1, T2> &pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

// Boost stdpair hasher from <boost/container_hash/hash.hpp>
// Copyright 2005-2014 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Based on Peter Dimov's proposal
//  http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1756.pdf
//  issue 6.18.
//
//  This also contains public domain code from MurmurHash. From the
//  MurmurHash header:

// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.
struct pair_hash_boost
{
private:
  template <typename SizeT>
  inline void hash_combine_impl(SizeT& seed, SizeT value)
  {
    seed ^= value + 0x9e3779b9 + (seed<<6) + (seed>>2);
  }
public:
  template <class T1, class T2>
  std::size_t operator() (const std::pair<T1, T2> &pair) const {
    std::size_t seed = 0;
    hash_combine_impl(seed, pair.first);
    hash_combine_impl(seed, pair.second);
    return seed;
  }


};
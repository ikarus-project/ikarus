// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once
// since python does not support passing python float by reference to a double&, we have to wrap everything
// see also https://pybind11.readthedocs.io/en/stable/faq.html#limitations-involving-reference-arguments
template <typename T>
struct ValueWrapper
{
  T val;
  ValueWrapper operator+(const ValueWrapper& v) const { return ValueWrapper{val + v.val}; }
  ValueWrapper operator-(const ValueWrapper& v) const { return ValueWrapper{val - v.val}; }
  ValueWrapper operator-() const { return ValueWrapper{-val}; }
  ValueWrapper operator*(T value) const { return ValueWrapper{val * value}; }
  ValueWrapper& operator+=(const ValueWrapper& v) {
    val += v.val;
    return *this;
  }
  ValueWrapper& operator*=(T v) {
    val *= v;
    return *this;
  }

  friend ValueWrapper operator*(T f, const ValueWrapper& v) { return ValueWrapper{f * v.val}; }
};
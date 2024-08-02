// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#define ENUM_BINDINGS(Type)                                        \
  py::enum_<Type> Type##Enum(m, #Type);                            \
  Type Type##EnumV = Type::BEGIN;                                  \
  Ikarus::increment(Type##EnumV);                                  \
  for (; Type##EnumV != Type::END; Ikarus::increment(Type##EnumV)) \
    Type##Enum.value(toString(Type##EnumV).c_str(), Type##EnumV);

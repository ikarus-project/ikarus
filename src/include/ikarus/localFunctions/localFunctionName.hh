/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */



#pragma once
#include <regex>
#include <string>

#include <dune/common/classname.hh>

namespace Ikarus {

  /** Pretty printing the name of a local function expression */
  template <typename LF>
  auto localFunctionName(const LF& lf) {
    std::string name = Dune::className(lf);

    std::regex regexp0("Ikarus::StandardLocalFunction<(([a-zA-Z0-9_:<, ]*>){10})");
    name = regex_replace(name, regexp0, "SLF");

    std::regex regexp1("Ikarus::ProjectionBasedLocalFunction<(([a-zA-Z0-9_:<, ]*>){10})");
    name = regex_replace(name, regexp1, "PBLF");

    std::regex regexp2("Ikarus::");
    name = regex_replace(name, regexp2, "");

    std::regex regexp3("LocalFunctionDot");
    name = regex_replace(name, regexp3, "Dot");

    std::regex regexp4("LocalFunctionNegate");
    name = regex_replace(name, regexp4, "Negate");

    std::regex regexp5("LocalFunctionScale");
    name = regex_replace(name, regexp5, "Scale");

    std::regex regexp6("LocalFunctionSum");
    name = regex_replace(name, regexp6, "Sum");

    std::regex regexp7("const");
    name = regex_replace(name, regexp7, "");

    std::regex regexp8("ConstantExpr");
    name = regex_replace(name, regexp8, "Constant");

    std::regex regexp9("SqrtExpr");
    name = regex_replace(name, regexp9, "Sqrt");

    std::regex regexp10("NormSquaredExpr");
    name = regex_replace(name, regexp10, "NormSquared");

    return name;
  }
}  // namespace Ikarus
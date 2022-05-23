//
// Created by Alex on 03.05.2022.
//

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

    return name;
  }
}  // namespace Ikarus
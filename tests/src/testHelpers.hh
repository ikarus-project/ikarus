//
// Created by Alex on 20.05.2021.
//

#pragma once

MATCHER_P2(EigenApproxEqual, expect, prec,
           std::string(negation ? "isn't" : "is") + " approx equal to\n" + ::testing::PrintToString(expect)
               + "\nwith precision " + ::testing::PrintToString(prec)) {
  if constexpr (requires { arg.isApprox(expect, prec); })
    return arg.isApprox(expect, prec);
  else  // Eigen::DiagonalMatrix branch
    return arg.diagonal().isApprox(expect.diagonal(), prec);
}

MATCHER_P(EigenExactEqual, expect,
          std::string(negation ? "isn't" : "is") + " equal to" + ::testing::PrintToString(expect)) {
  return ((arg == expect) == true).all();
}

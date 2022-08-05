

#pragma once

MATCHER_P2(EigenApproxEqual, expect, prec,
           std::string(negation ? "isn't" : "is") + " approx equal to\n" + ::testing::PrintToString(expect)
               + "\nwith precision " + ::testing::PrintToString(prec)) {
  if constexpr (requires {
                  arg.isApprox(expect, prec);
                  (arg - expect).isMuchSmallerThan(1, prec);
                })
    return arg.isApprox(expect, prec) or (arg - expect).isZero(prec);
  else if constexpr (requires { arg.isApprox(expect, prec); })
    return arg.isApprox(expect, prec);
  else  // Eigen::DiagonalMatrix branch
    return arg.diagonal().isApprox(expect.diagonal(), prec) or (arg.diagonal() - expect.diagonal()).isZero(prec);
}

MATCHER_P(EigenExactEqual, expect,
          std::string(negation ? "isn't" : "is") + " equal to" + ::testing::PrintToString(expect)) {
  return ((arg == expect) == true).all();
}

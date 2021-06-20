//
// Created by Alex on 20.05.2021.
//

#pragma once

MATCHER_P2(EigenApproxEqual, expect, prec,
           std::string(negation ? "isn't" : "is") + " approx equal to"
               + ::testing::PrintToString(expect) + "\nwith precision "
               + ::testing::PrintToString(prec)) {
  return arg.isApprox(expect, prec);
}

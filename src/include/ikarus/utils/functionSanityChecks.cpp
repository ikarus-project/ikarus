//
// Created by alex on 3/8/22.
//

#include "functionSanityChecks.hh"

#include "findLineSegment.hh"

#include <matplot/matplot.h>
#include <matplot/util/colors.h>

namespace Ikarus {
  double drawResultAndReturnSlope(std::string&& functionName, const std::function<double(double)>& ftfunc, bool draw) {
    using namespace matplot;
    std::vector<double> t = logspace(-8, 0, 100);
    Eigen::Map<Eigen::VectorXd> tE(t.data(), t.size());
    std::vector<double> ftevaluated = transform(t, ftfunc);
    Eigen::Map<Eigen::VectorXd> yE(ftevaluated.data(), ftevaluated.size());

    std::vector<double> fexpectedSlope = transform(t, [](auto t) { return t * t; });
    const int rangeSize                = 10;
    const auto [poly, range]           = Ikarus::findLineSegment(tE.array().log10(), yE.array().log10(), rangeSize);

    if (draw) {
      auto f   = figure(true);
      auto ax1 = gca();
      hold(ax1, true);
      std::vector<double> tOfRange(rangeSize);
      std::vector<double> fInRange(rangeSize);
      auto tET = tE(range);
      auto yET = yE(range);

      for (int i = 0; auto r : tET) {
        tOfRange[i] = r;
        fInRange[i] = yET[i];
        ++i;
      }

      auto l0          = ax1->loglog(t, ftevaluated);
      auto lexpected   = ax1->loglog(t, fexpectedSlope, "--");
      auto lFoundRange = ax1->loglog(tOfRange, fInRange);
      l0->line_width(2);
      lexpected->line_width(2);
      lFoundRange->line_width(4);
      lFoundRange->color("magenta");
      l0->color("blue");
      lexpected->color("red");
      xlabel("h");
      ylabel("Approximation error ");
      title(functionName + "check");
      f->show();
    }

    return poly.coefficients()[1];
  }

}  // namespace Ikarus
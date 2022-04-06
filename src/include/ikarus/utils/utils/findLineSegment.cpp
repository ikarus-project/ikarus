//
// Created by alex on 3/18/22.
//
#include "findLineSegment.hh"

namespace Ikarus {
  /*
   * This function is inspired from
   * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/identify_linear_piece.m
   */
  std::tuple<Dune::Functions::Polynomial<double>, decltype(Eigen::seq(0,0))> findLineSegment(const Eigen::VectorXd& x, const Eigen::VectorXd& y, int segmentSize) {
    Eigen::VectorXd errors = Eigen::VectorXd::Zero(x.size() - segmentSize);
    std::vector<Dune::Functions::Polynomial<double>> lines;
    for (int i = 0; i < errors.size(); ++i) {
      auto range = Eigen::seq(i, i + segmentSize);

      auto [poly, error] = polyfit(x(range), y(range), 1);
      errors(i)          = error;
      lines.push_back(poly);
    }
    auto minEle = std::ranges::min_element(errors.begin(), errors.end());
    auto index  = std::distance(errors.begin(), minEle);
    auto range  = Eigen::seq(index, index + segmentSize);
    return std::make_tuple(lines[index], range);
  }
}  // namespace Ikarus
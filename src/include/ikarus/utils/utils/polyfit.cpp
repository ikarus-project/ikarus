//
// Created by lex on 08/03/2022.
//
#include "polyfit.h"
namespace Ikarus {
  std::tuple<Dune::Functions::Polynomial<double>,double> polyfit(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>  &y, const int deg) {
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(x.size(), deg + 1);
    for (int j = 1; j < deg + 1; ++j)
      A.col(j) = A.col(j - 1).cwiseProduct(x);

    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(A.rows(), A.cols());
    qr.compute(A);
    Eigen::VectorXd coeffs = qr.solve(y);

    std::vector<double> coeffsSTD(coeffs.begin(), coeffs.end());
    Dune::Functions::Polynomial<double> poly(std::move(coeffsSTD));
    return std::make_tuple(poly,(A*coeffs-y).norm());
  }
}  // namespace Ikarus
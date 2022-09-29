//
#include <config.h>
#include <dune/common/test/testsuite.hh>
#include <dune/common/parallel/mpihelper.hh>
using Dune::TestSuite;
//#include "testHelpers.hh"
//
//#include <matplot/matplot.h>
//
//#include <Eigen/Core>
//
//#include <ikarus/utils/polyfit.hh>
//using namespace Catch;
//
//TEST_CASE("PolyFitTest: PolyFitTest1", "[polyFitTest.cpp]") {
//  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(10, 0, 10);
//  Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(10, 2, 20);
//
//  auto [poly, normE] = Ikarus::polyfit(x, y, 1);
//  CHECK(2.0 == Catch::Approx(poly.coefficients()[0]));
//  CHECK(1.8 == Catch::Approx(poly.coefficients()[1]));
//  CHECK(1e-14 > normE);
//}
//
//TEST_CASE("PolyFitTest: PolyFitTest2", "[polyFitTest.cpp]") {
//  const double factor = 7.6;
//  Eigen::VectorXd x   = Eigen::VectorXd::LinSpaced(10, 0, 10);
//  Eigen::VectorXd y   = 7 * x.array().cwiseProduct(x.array()).matrix();
//  for (int i = 0; i < y.size(); ++i) {
//    y[i] += (1 - i / 10.0) * factor - (1 - i * i / 10.0) * factor + std::sin(i / 10.0);
//  }
//
//  auto [poly, normE] = Ikarus::polyfit(x, y, 2);
//
//  CHECK(-0.0038062785674569739 == Catch::Approx(poly.coefficients()[0]));
//  CHECK(-0.58760441700969401 == Catch::Approx(poly.coefficients()[1]));
//  CHECK(7.6138682871655829 == Catch::Approx(poly.coefficients()[2]));
//  CHECK(0.0082367593944499204 == Catch::Approx(normE));
//}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    TestSuite t;

    //t.subTest(SimpleAssemblersTest());

    return t.exit();
}
//
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
using Dune::TestSuite;
//
// #include "testHelpers.hh"
//
// #include <array>
// #include <fstream>
// #include <unordered_set>
// #include <vector>
//
// #include <Eigen/Core>
//
// #include <ikarus/manifolds/realTuple.hh>
// #include <ikarus/manifolds/unitVector.hh>
// #include <ikarus/utils/algorithms.hh>
//
// static constexpr double tol = 1e-15;
//
// TEST_CASE("DefaultVariableTest: UnitVectorDirector", "[manifoldsTest.cpp]") {
//  using namespace Ikarus;
//  UnitVector<double, 3> a{UnitVector<double, 3>::CoordinateType::UnitZ()};
//  a.update(Eigen::Vector<double, 2>::UnitX());
//  const auto aExpected = Eigen::Vector<double, 3>(1.0 / sqrt(2), 0.0, 1.0 / sqrt(2));
//  CHECK_THAT(a.getValue(), EigenApproxEqual(aExpected, tol));
//
//  a = update(a, Eigen::Vector<double, 2>::UnitY());
//  CHECK_THAT(a.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(0.5, 1.0 / sqrt(2), 0.5), tol));
//
//  const auto d = update(a, Eigen::Vector<double, 2>::UnitY());
//
//  CHECK_THAT(d.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(0.18688672392660707344, 0.97140452079103167815,
//                                                                     -0.14644660940672621363),
//                                            tol));
//
//  UnitVector<double, 3> b{a};
//
//  CHECK_THAT(b.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(0.5, 1.0 / sqrt(2), 0.5), tol));
//
//  UnitVector<double, 3> c{UnitVector<double, 3>{Eigen::Vector<double, 3>::UnitZ() * 2.0}};  // move constructor test
//  CHECK_THAT(c.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(0.0, 0.0, 1.0), tol));
//
//  c.setValue(Eigen::Vector<double, 3>(13.0, -5.0, 1.0));
//  CHECK_THAT(c.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(13.0, -5.0, 1.0).normalized(), tol));
//
//  b = a;
//  CHECK(a == b);
//  const auto testVec = Eigen::Vector<double, 3>(127.0, -5.0, 1.0);
//  b.setValue(testVec);
//
//  CHECK_THAT(b.getValue(), EigenApproxEqual(testVec.normalized(), tol));
//
//  auto e{std::move(a)};
//
//  e.setValue(Eigen::Vector<double, 3>(0.0, 0.0, -1.0));
//
//  e = update(e, Eigen::Vector<double, 2>::UnitY());
//
//  const auto eExpected = Eigen::Vector<double, 3>(0, 1.0 / sqrt(2), -1.0 / sqrt(2));
//  CHECK_THAT(e.getValue(), EigenApproxEqual(eExpected, tol));
//}

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  // t.subTest(SimpleAssemblersTest());

  return t.exit();
}
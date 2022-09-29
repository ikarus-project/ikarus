//
//
#include <config.h>
#include <dune/common/test/testsuite.hh>
#include <dune/common/parallel/mpihelper.hh>
using Dune::TestSuite;
//
//#include "testHelpers.hh"
//
//#include <ikarus/utils/autodiffHelper.hh>
//
//template <typename Scalar>
//Eigen::Vector<Scalar, 2> f(const Eigen::Vector<Scalar, 3>& x) {
//  return (x.array() * (x.array().sin())).template segment<2>(0);
//}
//
//TEST_CASE("AutoDiffHelper: hessianN", "[testAutodiffHelper.cpp]") {
//  Eigen::Vector3d xd;
//  xd << 1.0, 2.0, 3.0;
//  Eigen::Vector3dual2nd x = xd;
//  Eigen::Vector2dual2nd u;
//  std::array<Eigen::Vector<double, 3>, 2> g;
//  std::array<Eigen::Matrix<double, 3, 3>, 2> h;
//  Ikarus::hessianN(f<autodiff::dual2nd>, wrt(x), at(x), u, g, h);
//
//  for (int i = 0; i < 2; ++i) {
//    Eigen::Vector3d gExpected;
//    Eigen::Matrix3d hExpected;
//    gExpected.setZero();
//    hExpected.setZero();
//    gExpected[i]    = sin(xd[i]) + cos(xd[i]) * xd[i];
//    hExpected(i, i) = 2 * +cos(xd[i]) - xd[i] * sin(xd[i]);
//    CHECK_THAT(g[i], EigenApproxEqual(gExpected, 1e-14));
//    CHECK_THAT(h[i], EigenApproxEqual(hExpected, 1e-14));
//  }
//}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    TestSuite t;

    //t.subTest(SimpleAssemblersTest());

    return t.exit();
}
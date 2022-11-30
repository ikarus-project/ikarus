//
//
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
using Dune::TestSuite;

#include <ikarus/utils/duneUtilities.hh>

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  autodiff::real valDual = 7.0;
  valDual[1]             = 1;
  Python::start();
  auto pyLambda = Python::Conversion<autodiff::Real<1, double>>::toPy(valDual);

  autodiff::real valExpected;
  Python::Conversion<autodiff::Real<1, double>>::toC(pyLambda, valExpected);

  t.check(valDual == valExpected);

  return t.exit();
}
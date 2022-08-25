
#include <config.h>

#include <catch2/catch_session.hpp>

#include <dune/common/parallel/mpihelper.hh>

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
  int result = Catch::Session().run(argc, argv);
  return result;
}
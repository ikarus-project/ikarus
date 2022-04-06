//
// Created by Alex on 21.04.2021.
//
#include <config.hh>

#include <gmock/gmock.h>

#include <dune/common/parallel/mpihelper.hh>

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);

  std::cout << "Hello from ikarus!" << std::endl;
}

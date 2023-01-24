# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

#message("Find MPI: ")
#find_package(MPI QUIET)
#message("Find METIS: ")
#find_package(METIS REQUIRED)
#include(AddMETISFlags)
#message("Find SuperLU: ")
#find_package(SuperLU REQUIRED)
#message("Find QuadMath: ")
#find_package(QuadMath)
#include(AddSuperLUFlags)



target_link_dune_default_libraries(${PROJECT_NAME}) # link compiled dune libs
add_dune_all_flags(${PROJECT_NAME})
add_dune_pythonlibs_flags(${PROJECT_NAME})

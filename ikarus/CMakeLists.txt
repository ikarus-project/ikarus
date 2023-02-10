# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.18)

file(
  GLOB_RECURSE sources
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS ../ikarus/*.cpp
)

dune_add_library(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} PUBLIC ${DUNE_LIBS})
target_link_dune_default_libraries(${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

add_subdirectory(assembler)
add_subdirectory(controlRoutines)
add_subdirectory(finiteElements)
add_subdirectory(linearAlgebra)
add_subdirectory(solver)
add_subdirectory(utils)
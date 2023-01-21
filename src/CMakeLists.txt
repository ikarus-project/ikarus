# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.18)

file(
  GLOB_RECURSE sources
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS include/ikarus/*.cpp
)

dune_add_library(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} PUBLIC ${DUNE_LIBS})
target_link_dune_default_libraries(${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
#target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

include(../cmake/AddDependencies.cmake)

# make header path available
target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include>
)

add_subdirectory(include)

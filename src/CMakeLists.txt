# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-2.1-or-later

cmake_minimum_required(VERSION 3.18)

file(
  GLOB_RECURSE sources
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS include/ikarus/*.cpp
)

project(
  ikarus
  VERSION 0.2
  LANGUAGES C CXX
)

dune_add_library(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} PUBLIC ${DUNE_LIBS})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

include(../cmake/AddDependencies.cmake)

# make header path available
target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include>
)

add_subdirectory(include)

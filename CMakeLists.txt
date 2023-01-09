# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(
  ikarus
  VERSION 0.3
  LANGUAGES C CXX
)

add_definitions(
  -Wall
  -Wextra
  -Wshadow
  -Wnon-virtual-dtor
  -pedantic
  -Wpedantic
  -Wimplicit-fallthrough
  -Wnull-dereference
  -Wcast-align
  -Wunused
  -Woverloaded-virtual
  -Wmisleading-indentation
  -Wduplicated-cond
  -Wduplicated-branches
  -m64
)
add_definitions(-DDUNE_LOCALFEFUNCTIONS_USE_EIGEN=1)

# find dune packages to make cmake modules available
find_package(dune-common REQUIRED EXCLUDE_FROM_ALL)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${dune-common_MODULE_PATH})
include(DuneMacros)

dune_project()
dune_enable_all_packages()

add_subdirectory(src)
add_subdirectory(sandbox)
add_subdirectory(tests)
add_subdirectory(cmake/FormatTarget)

if(BUILD_DOCS)
  message("Build docs locally target ENABLED")
  add_subdirectory(docs/BuildLocally)
endif()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

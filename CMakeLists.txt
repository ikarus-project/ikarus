cmake_minimum_required(VERSION 3.18 FATAL_ERROR)


project(
  Ikarus
  VERSION 0.1
  LANGUAGES C CXX
)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

add_definitions(-Wall -Wpedantic -Wextra -m64)

#find dune packages to make cmake modules available
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
        ${dune-common_MODULE_PATH})
find_package(dune-fufem REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
        ${dune-fufem_MODULE_PATH})
message("${dune-fufem_MODULE_PATH}")
message("${dune-common_MODULE_PATH}")
message("${CMAKE_MODULE_PATH}")
include(DuneMacros)
include(DuneFufemMacros)
dune_project()

dune_enable_all_packages()

find_package(SuiteSparse OPTIONAL_COMPONENTS CHOLMOD UMFPACK)
include(AddSuiteSparseFlags)


add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(benchmarks)
enable_testing()
add_subdirectory(tests)
add_subdirectory(cmake/FormatTarget)

if(BUILD_DOCS)
  message("Build docs locally target ENABLED")
  add_subdirectory(docs/BuildLocally)
endif()
# finalize the dune project, e.g. generating config.h etc.

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
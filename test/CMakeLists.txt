cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(IkarusTests LANGUAGES CXX)

# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)
CPMAddPackage("gh:TheLartians/Format.cmake@1.7.0")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("debug mode")
  list(APPEND CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_debug/)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message("release mode")
  list(APPEND CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_release/)
endif()
include_directories(../libAddons/eigen/eigenAddons)
if(TEST_INSTALLED_VERSION)
  find_package(Ikarus REQUIRED)
else()
  CPMAddPackage(NAME Ikarus SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(IkarusTests ${sources})

cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH MY_PARENT_DIR)
list(APPEND CMAKE_PREFIX_PATH ${MY_PARENT_DIR}/cmake/modules)
message("${CMAKE_MODULE_PATH}")
#if(MINGW)
find_package(GTest CONFIG REQUIRED)
#else()
#find_package(Gtest REQUIRED)
#endif()
target_link_libraries(IkarusTests PRIVATE GTest::gmock_main PRIVATE Ikarus::Ikarus)
set_target_properties(IkarusTests PROPERTIES CXX_STANDARD 20)

# enable compiler warnings
if(NOT TEST_INSTALLED_VERSION)
    target_compile_options(Ikarus PUBLIC -Wall -Wpedantic -Wextra -Werror)
endif()

enable_testing()

include(GoogleTest)
gtest_add_tests(IkarusTests SOURCES ${sources})

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
  target_compile_options(Ikarus PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
  target_link_options(Ikarus PUBLIC -fprofile-arcs -ftest-coverage)
endif()

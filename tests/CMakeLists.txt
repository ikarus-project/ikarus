cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(IkarusTests LANGUAGES CXX)
# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

include(../cmake/tools.cmake)

list(APPEND CMAKE_PREFIX_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_release/
)

include_directories(../libAddons/eigen/eigenAddons)

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
message("GTest")
find_package(GTest CONFIG REQUIRED)
message("GTest")
add_executable(IkarusTests ${sources})

target_link_libraries(
  IkarusTests
  PRIVATE GTest::gmock_main
  PRIVATE Ikarus::Ikarus
)
set_target_properties(IkarusTests PROPERTIES CXX_STANDARD 20)

target_compile_options(IkarusTests PUBLIC -Wall -Wpedantic -Wextra -Werror)
target_compile_features(IkarusTests PRIVATE cxx_std_20)

include(GoogleTest)
gtest_discover_tests(IkarusTests)

if(ENABLE_TEST_COVERAGE)
  target_compile_options(IkarusTests PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
  target_link_options(IkarusTests PUBLIC -fprofile-arcs -ftest-coverage)
endif()

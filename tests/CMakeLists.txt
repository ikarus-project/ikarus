cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(IkarusTests LANGUAGES C CXX)

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

set(IWYU_ARGS "-Xiwyu --mapping_file=iwyu.imp")
include(../cmake/tools.cmake)

if(MINGW)
list(APPEND CMAKE_PREFIX_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_release/
)
endif()

include_directories(../libAddons/eigen/eigenAddons)

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/LoadControlTestWithDuneGridAndBasis.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/AssemblerTest.cpp)

find_package(GTest CONFIG REQUIRED)
add_executable(IkarusTests ${sources})

target_link_libraries(
  IkarusTests
  PRIVATE GTest::gmock
  PUBLIC Ikarus::Ikarus
)
set_target_properties(IkarusTests PROPERTIES CXX_STANDARD 20)

target_compile_options(
  IkarusTests PUBLIC -Wall -Wpedantic -Wextra
                      -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-deprecated-copy)
target_compile_options(IkarusTests PUBLIC -fconcepts-diagnostics-depth=5)

if(MINGW)
  target_compile_options(IkarusTests PUBLIC -Wa,-mbig-obj -ftrack-macro-expansion=0)
endif()

target_compile_features(IkarusTests PRIVATE cxx_std_20)

include(GoogleTest)
gtest_discover_tests(IkarusTests)

if(ENABLE_TEST_COVERAGE)
  message("Enable Test Coverage")
  target_compile_options(IkarusTests PUBLIC -O0 -fprofile-arcs -ftest-coverage -fPIC)
  target_link_options(IkarusTests PUBLIC -fprofile-arcs -ftest-coverage)
endif()

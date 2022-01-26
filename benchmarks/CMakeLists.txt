cmake_minimum_required(VERSION 3.18)
message("Configure googlebenchmarks: ")
project(IkarusBenchmarks LANGUAGES C CXX)
file(GLOB benchSourceFiles CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_executable(IkarusBenchmarks ${benchSourceFiles})

list(APPEND CMAKE_PREFIX_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_release/
)
include_directories(../libAddons/eigen/eigenAddons)
find_package(benchmark REQUIRED)

target_link_libraries(
  IkarusBenchmarks
  PRIVATE Ikarus::Ikarus
  PRIVATE benchmark::benchmark_main
)
target_compile_features(IkarusBenchmarks PRIVATE cxx_std_20)

cmake_minimum_required(VERSION 3.16)
message("Configure googlebenchmarks: ")
project(IkarusBenchmarks LANGUAGES CXX)
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
file(GLOB benchSourceFiles CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_executable(IkarusBenchmarks ${benchSourceFiles})

list(APPEND CMAKE_PREFIX_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_release/
)

find_package(benchmark REQUIRED)

target_link_libraries(
  IkarusBenchmarks
  PRIVATE Ikarus::Ikarus
  PRIVATE benchmark::benchmark_main
)
target_compile_features(IkarusBenchmarks PRIVATE cxx_std_20)

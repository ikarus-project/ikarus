cmake_minimum_required(VERSION 3.18)

file(
  GLOB_RECURSE sources
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS include/ikarus/*.cpp
)

project(
  Ikarus
  VERSION 0.1
  LANGUAGES C CXX
)

dune_add_library(${PROJECT_NAME} ${sources})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

include(../cmake/AddDependencies.cmake)

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include>
)

add_library(Ikarus::Ikarus ALIAS Ikarus)
file(
  GLOB_RECURSE sources
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS include/ikarus/*.cpp
)
message("${sources}")
cmake_minimum_required(VERSION 3.18)
project(
  Ikarus
  VERSION 1.0
  LANGUAGES CXX
)

add_library(${PROJECT_NAME} ${sources})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

list(APPEND CMAKE_PREFIX_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_release/
)
include(../cmake/AddDependencies.cmake)

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC Eigen3::Eigen
  PUBLIC spdlog::spdlog
  PUBLIC dunecommon
  PUBLIC dunegeometry
  PUBLIC dunegrid
  PUBLIC ${SuiteSparse_LIBRARIES}
  # PUBLIC muesli
  PUBLIC Matplot++::matplot
  # PUBLIC Python3::Module PRIVATE pybind11::headers
)

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include>
)

add_library(Ikarus::Ikarus ALIAS Ikarus)

include(GNUInstallDirs)
install(TARGETS Ikarus EXPORT IkarusTargets)
# install_targets(Ikarus) install(DIRECTORY include/ikarus DESTINATION include COMPONENT Devel )

# configure_and_install(${CMAKE_SCRIPTS_ROOT}/simple-config.cmake.in Ikarus SameMajorVersion)

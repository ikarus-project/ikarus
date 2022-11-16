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
#
add_subdirectory(include)
#add_library(Ikarus::ikarus ALIAS ikarus)
#
#include(GNUInstallDirs)
#set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ikarus)
#
#install(
#  TARGETS ikarus
#  EXPORT IkarusTargets
#  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#)
#
#set_target_properties(ikarus PROPERTIES EXPORT_NAME Ikarus)
#
#install(
#  EXPORT IkarusTargets
#  FILE IkarusTargets.cmake
#  NAMESPACE Ikarus::
#  DESTINATION ${INSTALL_CONFIGDIR}
#)
#
#include(CMakePackageConfigHelpers)
#configure_package_config_file(
#        ${CMAKE_CURRENT_LIST_DIR}/../cmake/IkarusConfig.cmake.in
#        ${CMAKE_CURRENT_BINARY_DIR}/IkarusConfig.cmake INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
#)
#
#
#write_basic_package_version_file(
#  ${CMAKE_CURRENT_BINARY_DIR}/IkarusConfigVersion.cmake
#  VERSION ${PROJECT_VERSION}
#  COMPATIBILITY AnyNewerVersion
#)
#
#
## Install the config, configversion and custom find modules
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/IkarusConfig.cmake
#              ${CMAKE_CURRENT_BINARY_DIR}/IkarusConfigVersion.cmake
#        DESTINATION ${INSTALL_CONFIGDIR}
#)
#
## Exporting from the build tree
#export(
#  EXPORT IkarusTargets
#  FILE ${CMAKE_CURRENT_BINARY_DIR}/IkarusTargets.cmake
#  NAMESPACE Ikarus::
#)
#
## Register package in the User Package Registry
#export(PACKAGE Ikarus)

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
project(IkarusFormat)

if(ADD_FORMATTARGET)
  include(../CPM.cmake)
  CPMAddPackage(
    NAME Format.cmake
    VERSION 1.7.0
    GITHUB_REPOSITORY TheLartians/Format.cmake
    OPTIONS "FORMAT_SKIP_CMAKE NO" "CMAKE_FORMAT_EXCLUDE ../modules/FindSuiteSparse.cmake"
  )
endif()

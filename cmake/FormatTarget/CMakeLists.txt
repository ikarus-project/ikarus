cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(IkarusFormat)

if(ADD_FORMATTARGET)
  include(../CPM.cmake)
  CPMAddPackage(
    NAME Format.cmake
    VERSION 1.7.3
    GITHUB_REPOSITORY TheLartians/Format.cmake
    OPTIONS "FORMAT_SKIP_CMAKE NO"
  )
endif()

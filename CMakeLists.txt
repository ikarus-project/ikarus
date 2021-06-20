cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

# ---- Project ----

# Note: update this to your new project's name and version
project(
  Ikarus
  VERSION 1.0
  LANGUAGES CXX
)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# ---- Add dependencies via CPM ----
# see https://github.com/TheLartians/CPM.cmake for more info

include(cmake/CPM.cmake)

# PackageProject.cmake will be used to make our target installable
CPMAddPackage("gh:TheLartians/PackageProject.cmake@1.6.0")

# ---- Add source files ----

# Note: globbing sources is considered bad practice as CMake's generators may not detect new files
# automatically. Keep that in mind when changing files, or explicitly mention them here.
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp")

# ---- Create library ----
add_library(Ikarus ${headers} ${sources})

set_target_properties(Ikarus PROPERTIES CXX_STANDARD 20)

# being a cross-platform target, we enforce standards conformance on MSVC
target_compile_options(Ikarus PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("debug mode")
  list(APPEND CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_debug/)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message("release mode")
  list(APPEND CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/../../Ikarus_Dependencies/Dependencies_release/)
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

message("Find Eigen: ")
find_package(Eigen3 3.3.9 REQUIRED)
message("Find spdlog: ")
find_package(spdlog REQUIRED)
message("Find dune-common: ")
find_package(dune-common REQUIRED )
message("Find dune-geometry: ")
find_package(dune-geometry  REQUIRED)
message("Find dune-grid: ")
find_package(dune-grid REQUIRED)
message("Find SuiteSparse: ")
if(MINGW)
find_package(SuiteSparse CONFIG REQUIRED)
else()
find_package(SuiteSparse REQUIRED)
endif()

#message("Find muesli: ")
#find_package(muesli CONFIG REQUIRED)
#message("Find matplotc++: ")
#find_package(Matplot++ CONFIG REQUIRED)
message("Find PythonLibs: ")
find_package (Python3 COMPONENTS Interpreter Development)
message("${Python3_STDLIB}")
#message("Find pybind11: ")
#find_package(pybind11 CONFIG REQUIRED)


# Link dependencies
target_link_libraries(Ikarus  PUBLIC Eigen3::Eigen
        PUBLIC spdlog::spdlog
        PUBLIC dunecommon
        PUBLIC dunegeometry
        PUBLIC dunegrid
        PUBLIC ${SuiteSparse_LIBRARIES}
#        PUBLIC muesli
#        PUBLIC Matplot++::matplot
#        PUBLIC Python3::Module
#        PRIVATE pybind11::headers
)
target_include_directories(
  Ikarus PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                 $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)

packageProject(
  NAME ${PROJECT_NAME}
  VERSION ${PROJECT_VERSION}
  NAMESPACE ${PROJECT_NAME}
  BINARY_DIR ${PROJECT_BINARY_DIR}
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
  INCLUDE_DESTINATION include/${PROJECT_NAME}-${PROJECT_VERSION}
  VERSION_HEADER "${VERSION_HEADER_LOCATION}"
  COMPATIBILITY SameMajorVersion
)


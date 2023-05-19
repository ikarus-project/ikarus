# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
set(IWYU_ARGS "-Xiwyu --mapping_file=iwyu.imp")
include(../../cmake/tools.cmake)
find_dependency(Eigen3 3.3.9 REQUIRED)
find_dependency(spdlog REQUIRED)
find_dependency(autodiff REQUIRED)
find_dependency(Matplot++ REQUIRED)

file(GLOB programSourceFiles CONFIGURE_DEPENDS *.cpp)

foreach(programSourceFile ${programSourceFiles})
  get_filename_component(programName ${programSourceFile} NAME_WLE)
  dune_add_test(SOURCES ${programSourceFile} LINK_LIBRARIES ikarus LABELS cpp)
  target_compile_features(${programName} PUBLIC cxx_std_20)
  set_tests_properties(${programName} PROPERTIES TIMEOUT 1000)
  add_dune_pythonlibs_flags(${programName})
  find_package(PythonLibs REQUIRED)
  include_directories(${PYTHON_INCLUDE_DIRS})
  target_link_libraries(${programName} PUBLIC ${PYTHON_LIBRARIES})
  target_compile_options(${programName} PUBLIC -ftemplate-backtrace-limit=0)
  if(ENABLE_TEST_COVERAGE)
    message("Enable Test Coverage")
    target_compile_options(${programName} PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
    target_link_options(${programName} PUBLIC -fprofile-arcs -ftest-coverage)
  endif()

endforeach()

file(COPY testFiles/unstructuredTrianglesfine.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/testFiles)

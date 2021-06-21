cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(IkarusExamples LANGUAGES CXX)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)
CPMAddPackage(NAME Ikarus SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)

file(GLOB programSourceFiles CONFIGURE_DEPENDS source/*.cpp)
foreach(programSourceFile ${programSourceFiles})
  get_filename_component(programName ${programSourceFile} NAME_WLE)
  add_executable(${programName} ${programSourceFile})
  target_link_libraries(${programName} PRIVATE Ikarus::Ikarus)
  set_target_properties(${programName} PROPERTIES CXX_STANDARD 20 OUTPUT_NAME "Ikarus")
endforeach()

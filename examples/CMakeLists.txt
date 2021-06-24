cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(IkarusExamples LANGUAGES CXX)

file(GLOB programSourceFiles CONFIGURE_DEPENDS src/*.cpp)

foreach(programSourceFile ${programSourceFiles})
  get_filename_component(programName ${programSourceFile} NAME_WLE)
  add_executable(${programName} ${programSourceFile})
  target_link_libraries(${programName} PRIVATE Ikarus::Ikarus)
  set_target_properties(${programName} PROPERTIES CXX_STANDARD 20 OUTPUT_NAME "Ikarus")
endforeach()

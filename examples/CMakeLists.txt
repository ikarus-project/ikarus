cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(IkarusExamples LANGUAGES CXX)

file(GLOB programSourceFiles CONFIGURE_DEPENDS src/*.cpp)
include_directories(../libAddons/eigen/eigenAddons)
ADD_CUSTOM_TARGET(AllExamples)

foreach(programSourceFile ${programSourceFiles})
  get_filename_component(programName ${programSourceFile} NAME_WLE)
  add_executable(${programName} ${programSourceFile})
  target_link_libraries(${programName} PRIVATE Ikarus::Ikarus)
  set_target_properties(${programName} PROPERTIES CXX_STANDARD 20 OUTPUT_NAME "Ikarus")
  target_compile_options(
    ${programName} PUBLIC -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable
                          -Wno-deprecated-copy
  )
  ADD_DEPENDENCIES(AllExamples ${programName})
endforeach()

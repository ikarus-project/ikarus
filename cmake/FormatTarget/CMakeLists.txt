cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
project(IkarusFormat)


file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../../include/*.h")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../../source/*.cpp")

file(GLOB_RECURSE CMAKE_FILES
        CONFIGURE_DEPENDS CMakeLists.txt
        )

if(ADD_FORMATTARGET)
#    clang_format(Ikarus ${headers} ${sources})
#    file(GLOB_RECURSE CMAKE_FILES
#            CMakeLists.txt
#            )
#    cmake_format(cmake_format_Ikarus ${CMAKE_FILES})
    #  add_dependencies(format cmake_format_Ikarus)

        CPMAddPackage(
                NAME Format.cmake
                VERSION 1.7.0
                GITHUB_REPOSITORY TheLartians/Format.cmake
                OPTIONS
                "FORMAT_SKIP_CMAKE NO"
                "CMAKE_FORMAT_EXCLUDE cmake/"
        )
endif()




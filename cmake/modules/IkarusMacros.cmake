get_filename_component(IKARUS_CMAKE_DIR "IkarusMacros.cmake" PATH)
include(CMakeFindDependencyMacro)

find_dependency(spdlog)
find_dependency(Matplot++)
find_dependency(Eigen3)

set(IKARUS_lIBRARIES ikarus)

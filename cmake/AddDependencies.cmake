
message("====================")
if(MINGW)
  # cmake-format: off
  include(FetchContent)
  set(ikarusDepInstallDir ${CMAKE_SOURCE_DIR}/../Ikarus_Dependencies)
  FetchContent_Declare(
    ikarusDependencies
    URL https://github.com/IkarusRepo/IkarusDependencies/releases/download/v0.6.1/Dependencies_release.7z
    PREFIX        ${ikarusDepInstallDir}
    DOWNLOAD_DIR  ${ikarusDepInstallDir}/src
    SOURCE_DIR    ${ikarusDepInstallDir}/Dependencies_release
    INSTALL_COMMAND ""
    BUILD_COMMAND ""
    CONFIGURE_COMMAND ""
  )
  # cmake-format: on
  FetchContent_GetProperties(ikarusDependencies)
  if(NOT ikarusDependencies_POPULATED)
    FetchContent_Populate(ikarusDependencies)
  endif()
  FetchContent_MakeAvailable(ikarusDependencies)
endif()

message("Find MPI: ")
find_package(MPI QUIET)
#list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${dune-common_MODULE_PATH})
#include(DuneMacros)

message("Find METIS: ")
find_package(METIS REQUIRED)
include(AddMETISFlags)
message("Find SuperLU: ")
find_package(SuperLU REQUIRED)
include(AddSuperLUFlags)

add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
message("Find Eigen: ")
find_package(Eigen3 3.3.9 REQUIRED)
message("Find spdlog: ")
find_package(spdlog REQUIRED)
#message("Find dune-alugrid: ")
find_package(dune-alugrid REQUIRED)
#message("Find dune-foamgrid: ")
#find_package(dune-foamgrid REQUIRED)
#message("====================")
#message(${dune-alugrid_INCLUDE_DIRS})
#message("====================")
target_include_directories(${PROJECT_NAME} PUBLIC ${dune-alugrid_INCLUDE_DIRS})
message("====================")
message("${dune-alugrid_INCLUDE_DIRS}")
message("====================")
#
#message("Find dune-typetree: ")
#find_package(dune-typetree REQUIRED)
#message("Find dune-geometry: ")
#find_package(dune-geometry REQUIRED)
#message("Find dune-functions: ")
#find_package(dune-functions REQUIRED)
#message("Find dune-iga: ")
#find_package(dune-iga REQUIRED)
## message("Find dune-uggrid: ") find_package(dune-uggrid REQUIRED)
#message("Find dune-grid: ")
#find_package(dune-grid REQUIRED)
#message("Find dune-fufem: ")
#find_package(dune-fufem REQUIRED)
# message("Find muesli: ") find_package(muesli REQUIRED)


message("Find autodiff: ")
find_package(autodiff REQUIRED)

set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include(AddBLASLapackFlags)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
#find_package(dune-fufem REQUIRED)
#list(APPEND CMAKE_MODULE_PATH  ${dune-fufem_MODULE_PATH})
#message("${CMAKE_MODULE_PATH}")

include(AddPythonLibsFlags)
#message("Find SuiteSparse: ")
#if(MINGW OR MSVC)
#  find_package(SuiteSparse CONFIG REQUIRED CHOLMOD UMFPACK METIS)
#  set(SuiteSparseIncludeDirective ${SuiteSparse_LIBRARIES})
#else()
#  find_package(SuiteSparse OPTIONAL_COMPONENTS UMFPACK)
#  include(AddSuiteSparseFlags)
#
#  set(SuiteSparseIncludeDirective "SuiteSparse::SuiteSparse")
#endif()
message("Find matplotc++: ")
find_package(Matplot++ REQUIRED)
message("Find PythonLibs: ")




target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC Eigen3::Eigen
  PUBLIC SuperLU::SuperLU
  PUBLIC spdlog::spdlog
#  PUBLIC dunecommon
#  PUBLIC dunegeometry
  PUBLIC dunealugrid
#  PUBLIC dunegrid
#  PUBLIC dunefufem
#  PUBLIC duneuggrid
#  PUBLIC ${SuiteSparse_LIBRARIES}
  # PUBLIC muesli
  PUBLIC Matplot++::matplot
  PUBLIC autodiff::autodiff
  PUBLIC gfortran
)

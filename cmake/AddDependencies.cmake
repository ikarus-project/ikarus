# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
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
message("====================")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
  string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
          ${PROJECT_BINARY_DIR})
endif()
message("Find dune-common: ")
find_package(dune-common REQUIRED)

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
        ${dune-common_MODULE_PATH})

include(DuneMacros)

message("Find METIS: ")
find_package(METIS REQUIRED)
message("Find SuperLU: ")
find_package(SuperLU REQUIRED)
add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
message("Find Eigen: ")
find_package(Eigen3 3.3.9 REQUIRED)
message("Find spdlog: ")
find_package(spdlog REQUIRED)
message("Find dune-alugrid: ")
find_package(dune-alugrid REQUIRED)
message("Find dune-foamgrid: ")
find_package(dune-foamgrid REQUIRED)
message("====================")
message(${dune-alugrid_INCLUDE_DIRS})
message("====================")
target_include_directories(${PROJECT_NAME} PUBLIC ${dune-alugrid_INCLUDE_DIRS})



message("Find dune-typetree: ")
find_package(dune-typetree REQUIRED)
message("Find dune-geometry: ")
find_package(dune-geometry REQUIRED)
message("Find dune-functions: ")
find_package(dune-functions REQUIRED)
message("Find dune-iga: ")
find_package(dune-iga REQUIRED)
# message("Find dune-uggrid: ") find_package(dune-uggrid REQUIRED)
message("Find dune-grid: ")
find_package(dune-grid REQUIRED)

# message("Find muesli: ") find_package(muesli REQUIRED)

message("Find autodiff: ")
find_package(autodiff REQUIRED)
set(BLA_VENDOR OpenBLAS)
message("Find blas: ")
find_package(BLAS REQUIRED)
message("Find lapack: ")
find_package(LAPACK REQUIRED)
message("Find SuiteSparse: ")
if(MINGW OR MSVC)
  find_package(SuiteSparse CONFIG REQUIRED CHOLMOD UMFPACK METIS)
  set(SuiteSparseIncludeDirective ${SuiteSparse_LIBRARIES})
else()
  find_package(SuiteSparse REQUIRED CHOLMOD UMFPACK)
  set(SuiteSparseIncludeDirective "SuiteSparse::SuiteSparse")
endif()
message("Find matplotc++: ")
find_package(Matplot++ REQUIRED)
message("Find PythonLibs: ")
find_package(Python3 COMPONENTS Interpreter Development)
message("${Python3_STDLIB}")

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC Eigen3::Eigen
  PUBLIC METIS::METIS
  PUBLIC SuperLU::SuperLU
  PUBLIC LAPACK::LAPACK
  PUBLIC BLAS::BLAS
  PUBLIC spdlog::spdlog
  PUBLIC dunecommon
  PUBLIC dunegeometry
  PUBLIC dunealugrid
  PUBLIC dunegrid
  # PUBLIC duneuggrid
  PUBLIC ${SuiteSparseIncludeDirective}
  # PUBLIC muesli
  PUBLIC Matplot++::matplot
  PUBLIC autodiff::autodiff
  PUBLIC gfortran
)

get_target_property(OUT ${PROJECT_NAME} HEADER_FILE_ONLY)
message("====================")
message(STATUS ${OUT})
message(STATUS ${PROJECT_NAME})
message("====================")

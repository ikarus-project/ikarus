list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")

if(MINGW)
  # cmake-format: off
  include(FetchContent)
  set(ikarusDepInstallDir ${CMAKE_SOURCE_DIR}/../Ikarus_Dependencies)
  FetchContent_Declare(
    ikarusDependencies
    URL https://github.com/IkarusRepo/IkarusDependencies/releases/download/v0.41/Dependencies_release.7z
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
#list(APPEND CMAKE_PREFIX_PATH "/home/alex/Documents/dune/dune-alugrid/build-cmake/")
message("Find MPI: ")
find_package(MPI QUIET)

message("Find METIS: ")
find_package(METIS REQUIRED)
message("Find Eigen: ")
find_package(Eigen3 3.3.9 REQUIRED)
message("Find spdlog: ")
find_package(spdlog REQUIRED)
#message("Find alugrid: ")

list(APPEND CMAKE_PREFIX_PATH "/home/lex/Dokumente/dune/")

find_package(dune-alugrid REQUIRED)
message("====================")
message(${dune-alugrid_INCLUDE_DIRS})
message("====================")
target_include_directories(${PROJECT_NAME} PUBLIC ${dune-alugrid_INCLUDE_DIRS})
message("Find dune-functions: ")
find_package(dune-functions REQUIRED)
message("Find dune-common: ")
find_package(dune-common REQUIRED)
message("Find dune-geometry: ")
find_package(dune-geometry REQUIRED)
set(HAVE_UG true)
message("Find dune-uggrid: ")
find_package(dune-uggrid REQUIRED)
#message("${HAVE_UG}")
message("Find dune-grid: ")
find_package(dune-grid REQUIRED)


message("${HAVE_UG}")
message("${HAVE_UG}")
# message("Find muesli: ") find_package(muesli REQUIRED)

message("Find autodiff: ")
find_package(autodiff REQUIRED)

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
  PUBLIC spdlog::spdlog
  PUBLIC dunealugrid
  PUBLIC dunecommon
  PUBLIC dunegeometry
  PUBLIC dunegrid
  PUBLIC duneuggrid
  PUBLIC ${SuiteSparseIncludeDirective}
  # PUBLIC muesli
  PUBLIC Matplot++::matplot
  PUBLIC autodiff::autodiff
  PUBLIC gfortran
)
message("${HAVE_UG}")
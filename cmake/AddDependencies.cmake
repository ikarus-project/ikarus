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

message("Find MPI: ")
find_package(MPI QUIET)
message("Find Eigen: ")
find_package(Eigen3 3.3.9 REQUIRED)
message("Find spdlog: ")
find_package(spdlog REQUIRED)
message("Find dune-common: ")
find_package(dune-common REQUIRED)
message("Find dune-geometry: ")
find_package(dune-geometry REQUIRED)
message("Find dune-grid: ")
find_package(dune-grid REQUIRED)

#message("Find muesli: ")
#find_package(muesli REQUIRED)

message("Find autodiff: ")
find_package(autodiff REQUIRED)

message("Find SuiteSparse: ")
if(MINGW OR MSVC)
  find_package(SuiteSparse CONFIG REQUIRED CHOLMOD UMFPACK)
else()
  find_package(SuiteSparse REQUIRED CHOLMOD UMFPACK)
endif()
message("Find matplotc++: ")
find_package(Matplot++ REQUIRED)
message("Find PythonLibs: ")
find_package(Python3 COMPONENTS Interpreter Development)
message("${Python3_STDLIB}")

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC Eigen3::Eigen
  PUBLIC spdlog::spdlog
  PUBLIC dunecommon
  PUBLIC dunegeometry
  PUBLIC dunegrid
  PUBLIC ${SuiteSparse_LIBRARIES}
  PUBLIC UMFPACK
  PUBLIC CHOLMOD
  # PUBLIC muesli
  PUBLIC Matplot++::matplot
  PUBLIC autodiff::autodiff
  PUBLIC gfortran
)

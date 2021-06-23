list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")

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
message("Find SuiteSparse: ")
if(MINGW)
  find_package(SuiteSparse CONFIG REQUIRED)
else()
  find_package(SuiteSparse REQUIRED)
endif()
message("Find matplotc++: ")
find_package(Matplot++ REQUIRED)
message("Find PythonLibs: ")
find_package(Python3 COMPONENTS Interpreter Development)
message("${Python3_STDLIB}")

target_link_libraries(
  Ikarus
  PUBLIC Eigen3::Eigen
  PUBLIC spdlog::spdlog
  PUBLIC dunecommon
  PUBLIC dunegeometry
  PUBLIC dunegrid
  PUBLIC ${SuiteSparse_LIBRARIES}
  # PUBLIC muesli
  PUBLIC Matplot++::matplot
)

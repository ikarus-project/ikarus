message("Find MPI: ")
find_package(MPI QUIET)
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
message("Find autodiff: ")
find_package(autodiff REQUIRED)

message("Find matplotc++: ")
find_package(Matplot++ REQUIRED)
dune_register_package_flags(
        INCLUDE_DIRS
        ${Eigen3_INCLUDE_DIRS}
        ${spdlog_INCLUDE_DIRS}
        ${matplot_INCLUDE_DIRS}
        ${autodiff_INCLUDE_DIRS}

        LIBRARIES
        Eigen3::Eigen
        spdlog::spdlog
        Matplot++::matplot
        autodiff::autodiff
        )

target_link_dune_default_libraries(${PROJECT_NAME}) # link compiled dune libs
add_dune_all_flags(${PROJECT_NAME})
add_dune_pythonlibs_flags(${PROJECT_NAME})

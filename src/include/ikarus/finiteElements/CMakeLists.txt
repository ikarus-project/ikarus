# install headers
install(FILES feRequirements.hh feTraits.hh fEparameter.hh physicsHelper.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteElements
)

add_subdirectory(feBases)
add_subdirectory(mechanics)

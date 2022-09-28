#install headers
install(FILES
        algorithms.hh
        autodiffHelper.hh
        concepts.hh
        eigenDuneTransformations.hh
        findLineSegment.hh
        functionSanityChecks.hh
        linearAlgebraHelper.hh
        multiIndex.hh
        namedType.hh
        polyfit.hh
        traits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/utils)

add_subdirectory(drawing)
add_subdirectory(observer)
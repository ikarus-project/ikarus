set(IKARUS_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ikarus)

install(EXPORT ikarusTargets
        FILE ikarusTargets.cmake
        NAMESPACE ikarus::
        DESTINATION ${IKARUS_INSTALL_CONFIGDIR}
        COMPONENT cmake)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/ikarusConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
        ARCH_INDEPENDENT)

configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/IkarusConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/ikarusConfig.cmake
        INSTALL_DESTINATION ${IKARUS_INSTALL_CONFIGDIR}
        PATH_VARS IKARUS_INSTALL_CONFIGDIR)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/ikarusConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/ikarusConfigVersion.cmake
        DESTINATION ${IKARUS_INSTALL_CONFIGDIR})

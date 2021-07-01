add_custom_target(site
        COMMAND mkdocs build
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/site ${CMAKE_CURRENT_BINARY_DIR}/site
        COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_CURRENT_SOURCE_DIR}/site
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
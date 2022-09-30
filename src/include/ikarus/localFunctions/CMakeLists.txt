# install headers
install(
  FILES expressions.hh
        leafNodeCollection.hh
        localFunctionArguments.hh
        localFunctionHelper.hh
        localFunctionInterface.hh
        localFunctionName.hh
        meta.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/localFunctions
)

add_subdirectory(expressions)
add_subdirectory(impl)

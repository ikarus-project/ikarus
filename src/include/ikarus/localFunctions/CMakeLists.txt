# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-2.1-or-later

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

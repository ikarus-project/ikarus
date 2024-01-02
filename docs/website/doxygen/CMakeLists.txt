# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# mueller@ibb.uni-stuttgart.de SPDX-License-Identifier: LGPL-3.0-or-later
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS ../../ikarus/**/*.hh)
list(APPEND HEADERS modules.txt mainpage.txt)
message(STATUS ${HEADERS})
add_doxygen_target(DEPENDS ${HEADERS})
message(STATUS "Add doxygen target")

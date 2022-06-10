# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
add_subdirectory(ikarus)
set(CMAKE_CXX_STANDARD 20)

message(STATUS "PYPI_PACKAGE_NAME ${PYPI_PACKAGE_NAME}")
dune_python_configure_bindings(
  PATH
  "."
  PACKAGENAME
  ikarus
  CMAKE_METADATA_FLAGS
  DUNE_OPTS_FILE
  CMAKE_CXX_STANDARD
)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()

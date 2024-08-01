# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_python_targets(io __init__)

dune_add_pybind11_module(NAME _io SOURCES _io.cc)

set_property(
  TARGET _io
  PROPERTY LINK_LIBRARIES dunecommon ikarus
  APPEND
)

if(SKBUILD)
  install(TARGETS _io LIBRARY DESTINATION python/ikarus/io)
endif()

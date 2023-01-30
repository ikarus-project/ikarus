# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_python_targets(ikarus
  __init__
# toVoigt
)
#find_package(autodiff REQUIRED)
dune_add_pybind11_module(NAME _ikarus)
set_property(TARGET _ikarus PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid ikarus APPEND)
target_compile_options(_ikarus PRIVATE -fPIC)
target_compile_features(_ikarus PRIVATE cxx_std_20)
if(SKBUILD)
  install(TARGETS _ikarus LIBRARY DESTINATION python/ikarus)
endif()

# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES loadcontrol.hh
        loadcontrol.inl
        pathfollowing.hh
        pathfollowing.inl
        adaptivestepsizing.hh
        controlinfos.hh
        pathfollowingfunctions.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/controlroutines
)

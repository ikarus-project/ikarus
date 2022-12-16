# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.18)
project(Ikarus-Docs)
add_custom_target(
  site
  COMMAND xvfb-run -a mkdocs build --config-file mkdocs.insiders.yml
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)

cmake_minimum_required(VERSION 3.18)
project(Ikarus-Docs)
add_custom_target(
  localSite
  COMMAND mkdocs serve
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

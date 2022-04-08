cmake_minimum_required(VERSION 3.18)
project(Ikarus-Docs)
find_program(
  MKDOCS_EXECUTABLE
  NAMES mkdocs
  DOC "MkDocs documentation generation tool (http://www.mkdocs.org)" REQUIRED
)

add_custom_target(
  localSite
  COMMAND ${MKDOCS_EXECUTABLE} serve
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

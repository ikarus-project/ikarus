# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

configure_file(setpath.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setpath.py)

dune_python_add_test(
  NAME
  pylinearElasticTest
  SCRIPT
  linearElasticTest.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)

dune_python_add_test(
  NAME
  pynonLinearElasticTest
  SCRIPT
  nonLinearElasticTest.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)

# SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

configure_file(setpath.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setpath.py)

dune_python_add_test(
  NAME
  pylinearElasticTest
  SCRIPT
  linearelastictest.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)

dune_python_add_test(
  NAME
  pynonLinearElasticTest
  SCRIPT
  nonlinearelastictest.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)

dune_python_add_test(
  NAME
  pymaterial
  SCRIPT
  testmaterials.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)

if(HAVE_DUNE_IGA)
  dune_python_add_test(
    NAME
    pykltest
    SCRIPT
    kltest.py
    WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
    LABELS
    python
  )
endif()

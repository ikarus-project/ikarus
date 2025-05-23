# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

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
  pytrussTest
  SCRIPT
  trusstest.py
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

dune_python_add_test(
  NAME
  pydirichletvalues
  SCRIPT
  dirichletvaluetest.py
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}
  LABELS
  python
)

dune_python_add_test(
  NAME
  pyvtkwritertest
  SCRIPT
  vtkwritertest.py
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

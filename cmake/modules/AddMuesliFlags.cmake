# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_MUESLI for config.h
set(HAVE_MUESLI MUESLI)

# register all MUESLI related flags
if(HAVE_MUESLI)
  message("Muesli found yuhu")
  dune_register_package_flags(LIBRARIES  ${MUESLI} COMPILE_DEFINITIONS "ENABLE_MUESLI=1")
  else()
  message("Muesli not found crycry")
endif()



# add function to link against the MUESLI library
function(add_dune_MUESLI_flags _targets)
  if(HAVE_MUESLI)
    message("Muesli added")
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC  ${MUESLI})
      target_compile_definitions(${_target} PUBLIC ENABLE_MUESLI=1)

    endforeach(_target)
  endif()
endfunction(add_dune_MUESLI_flags)

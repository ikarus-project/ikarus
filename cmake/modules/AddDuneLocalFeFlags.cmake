# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# set HAVE_dune-localfefunctions for config.h
set(HAVE_DUNE_LOCALFEFUNCTIONS ${dune-localfefunctions_FOUND})

# register all dune-localfefunctions related flags
if(dune-localfefunctions_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_DUNE_LOCALFEFUNCTIONS=1;DUNE_LOCALFEFUNCTIONS_USE_EIGEN=1")
endif()

# add function to link against the dune-localfefunctions library
function(add_dune_dune-localfefunctions_flags _targets)
  if(dune-localfefunctions_FOUND)
    foreach(_target ${_targets})
      target_compile_definitions(${_target} PUBLIC ENABLE_DUNE_LOCALFEFUNCTIONS=1 PUBLIC DUNE_LOCALFEFUNCTIONS_USE_EIGEN=1)
    endforeach(_target)
  endif()
endfunction(add_dune_dune-localfefunctions_flags)

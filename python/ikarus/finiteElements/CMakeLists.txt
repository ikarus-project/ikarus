add_python_targets(ikarusa
  __init__
  _referenceelements
)
dune_add_pybind11_module(NAME _ikarusa)
set_property(TARGET _ikarusa PROPERTY LINK_LIBRARIES dunecommon dunegeometry APPEND)
if(SKBUILD)
    install(TARGETS _ikarusa LIBRARY DESTINATION python/ikarus/feRequirements)
endif()


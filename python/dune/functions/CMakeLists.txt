add_python_targets(ikarus
  __init__
)

find_package(autodiff REQUIRED)
dune_add_pybind11_module(NAME _ikarus)

set_property(TARGET _ikarus PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid APPEND)
if(SKBUILD)
    install(TARGETS _ikarus LIBRARY DESTINATION python/ikarus)
endif()

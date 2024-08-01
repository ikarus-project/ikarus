#include <config.h>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/io/vtkdatatag.hh>

//clang-format off
#define ENUM_BINDINGS(Type)                        \
  py::enum_<Type> Type##Enum(m, #Type);           \
  Type Type##EnumV = Type::BEGIN;                          \
  Ikarus::increment(Type##EnumV);                          \
  for (; Type##EnumV != Type::END; Ikarus::increment(Type##EnumV)) \
    Type##Enum.value(toString(Type##EnumV).c_str(), Type##EnumV);
// clang-format on

PYBIND11_MODULE(_io, m) {
  namespace py = pybind11;
  using namespace pybind11::literals;

  using namespace Ikarus::Vtk;
  using namespace Ikarus;
  ENUM_BINDINGS(DataTag);
}
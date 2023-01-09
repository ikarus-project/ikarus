#ifndef DUNE_PYTHON_FUNCTIONS_INTERPOLATE_HH
#define DUNE_PYTHON_FUNCTIONS_INTERPOLATE_HH

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>

#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune {
namespace Python {
namespace Functions {

template<class Basis, typename T>
void interpolate(const Basis& basis, pybind11::array_t<T> x, const std::function<T(typename Basis::GridView::template Codim<0>::Geometry::GlobalCoordinate)>& f)
{
  if (x.shape(0) != basis.size())
    throw std::runtime_error("Coefficient vector has wrong size");

  auto x1 = x.template mutable_unchecked<1>();

  using V = decltype(x1);
  auto x2 = Dune::Functions::HierarchicVectorWrapper<V, T>(x1);

  interpolate(basis, x2, f);
}

} /* namespace Functions */
} /* namespace Python */
} /* namespace Dune */

#endif

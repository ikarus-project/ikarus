import numpy as np

import dune.grid
from dune.grid import cartesianDomain
from dune.grid import yaspGrid
import dune.functions
from dune.generator import algorithm
import io

# we try to pass different dune-functions functions from the python side to C++
code="""
#include <dune/common/classname.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/python/functions/hierarchicvectorwrapper.hh>
#include <dune/python/functions/globalbasis.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

static const int dim = 2;
using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
using GV = typename Grid::LeafGridView;
using Signature = double(Dune::FieldVector<double,2>);
using GridViewFunction = Dune::Functions::GridViewFunction<Signature,GV>;

template<typename B, typename V, typename NTRE, typename R>
std::string callScalar(const Dune::Functions::DiscreteGlobalBasisFunction<B,V,NTRE,R> & f)
{
    return "DiscreteGlobalBasisFunction<...>";
}

std::string callScalar(const GridViewFunction & f)
{
    return "GridViewFunction<GV>";
}

std::string callVector(const std::vector<GridViewFunction> & f)
{
    return "vector<GridViewFunction<GV>>";
}
"""

class CppTypeInfo:
    def __init__(self,name,includes):
        self.cppTypeName = name
        self.cppIncludes = includes

scalarTypeName = "GridViewFunction"
vectorTypeName = "std::vector<GridViewFunction>"
callScalar = algorithm.load('callScalar', io.StringIO(code), CppTypeInfo(scalarTypeName,[]));
callVector = algorithm.load('callVector', io.StringIO(code), CppTypeInfo(vectorTypeName,[]));

# number of grid elements (in one direction)
gridSize = 4

# create a YaspGrid of the unit square
domain = cartesianDomain([0,0],[1,1],[gridSize,gridSize])
grid = yaspGrid(domain, dimgrid=2)

# create a nodal Lagrange FE basis of order 1 and order 2
basis1 = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=1))
basis2 = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=2))

# create a DOF vectors
N1 = len(basis1)
x1 = np.ndarray(N1)

N2 = len(basis2)
x2 = np.ndarray(N2)

# interpolate data
basis1.interpolate(x1, lambda x : np.linalg.norm(x-np.array([0.5,0.5]))-0.3)
basis2.interpolate(x2, lambda x : np.linalg.norm(x-np.array([0.5,0.5]))-0.3)

# create a grid function
f1 = basis1.asFunction(x1)
f2 = basis2.asFunction(x2)

# calling the concrete interface requires the lagrange basis header
incLagrange = "#include <dune/functions/functionspacebases/lagrangebasis.hh>\n"
callConcrete1 = algorithm.load('callScalar', io.StringIO(incLagrange + code), f1);
callConcrete2 = algorithm.load('callScalar', io.StringIO(incLagrange + code), f2);

print ("Calling into C++")
print ("Called as " + callConcrete1(f1))
print ("Called as " + callConcrete2(f2))
# try to call as GridViewFunction
print ("Called as " + callScalar(f1))
print ("Called as " + callScalar(f2))
print ("Called as " + callVector([f1,f1]))
print ("Called as " + callVector([f1,f2]))
print ("Called as " + callVector([f2,f2]))

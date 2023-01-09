from __future__ import absolute_import, division, print_function, unicode_literals

from .tree import Composite, DG, Lagrange, Power, Tree

duneFunctionsLayouts = {"lexicographic": "Lexicographic", "interleaved": "Interleaved"}

def indexMergingStrategy(blocked, layout):
    return "Dune::Functions::BasisBuilder::" + ("Blocked" if blocked else "Flat") + duneFunctionsLayouts[layout]


def preBasisTypeName(tree, gridViewTypeName):
    assert isinstance(tree, Tree)
    if isinstance(tree, Lagrange):
        scalarPreBasis = "Dune::Functions::LagrangePreBasis< " + gridViewTypeName + " , " + str(tree.order) + " >"
        if tree.dimRange != 1:
            IMS = indexMergingStrategy(False, "interleaved")
            return "Dune::Functions::PowerPreBasis< " + IMS + " , " + scalarPreBasis + " , " + str(tree.dimRange) + " >"
        else:
            return scalarPreBasis
    elif isinstance(tree, DG):
        raise Exception(repr(tree) + " not supported by dune-functions.")
    elif isinstance(tree, Composite):
        IMS = indexMergingStrategy(tree.blocked, tree.layout)
        ChildPreBases = " , ".join(preBasisTypeName(c, gridViewTypeName) for c in tree.children)
        return "Dune::Functions::CompositePreBasis< " + IMS + " , " + ChildPreBases + " >"
    elif isinstance(tree, Power):
        IMS = indexMergingStrategy(tree.blocked, tree.layout)
        ChildPreBasis = preBasisTypeName(tree.children[0], gridViewTypeName)
        return "Dune::Functions::PowerPreBasis< " + IMS + " , " + ChildPreBasis + " , " + str(tree.exponent) + " >"
    else:
        raise Exception("Unknown type of tree: " + repr(tree))


def defaultGlobalBasis(gridView, tree):
    from dune.functions import load

    headers = ["powerbasis", "compositebasis", "lagrangebasis", "subspacebasis", "defaultglobalbasis"]

    includes = []
    includes += list(gridView.cppIncludes)
    includes += ["dune/functions/functionspacebases/" + h + ".hh" for h in headers]

    typeName = "Dune::Functions::DefaultGlobalBasis< " + preBasisTypeName(tree, gridView.cppTypeName) + " >"

    return load(includes, typeName).GlobalBasis(gridView)

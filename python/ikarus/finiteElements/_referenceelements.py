from ..generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

def module():
    typeName = "Ikarus::FeRequirements<Dune::Geo::ReferenceElementImplementation<double," + str(dim) + "> >"
    includes = ["ikarus/python/geometry/referenceelements.hh"]
    typeHash = "referenceelements_" + hashIt(typeName)
    generator = SimpleGenerator("ReferenceElements", "Ikarus::Python")
    m = generator.load(includes, typeName, typeHash)
    return m

_duneReferenceElements = {}
def referenceElement(geometryType):
    try:
        geometryType = geometryType.type
    except:
        pass
    try:
        ref = _duneReferenceElements[geometryType]
    except KeyError:
        ref = module(geometryType.dim).general(geometryType)
        _duneReferenceElements[geometryType] = ref
    return ref

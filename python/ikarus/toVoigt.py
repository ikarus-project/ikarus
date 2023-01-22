
def  to_voigt() :
    from dune.common.hashit import hashIt
    typeName = "Ikarus::toVoigt< double ," + str(len(values)) + " >"
generator = SimpleGenerator("toVoigt", "Ikarus::Python")
    includes = includes + ["python/ikarus/toVoigt.hh"]
    typeHash = "tovoigt_" + hashIt(typeName)
    return generator.load(
        includes, typeName, typeHash,
        constructors, methods, bufferProtocol=True
    )


    def _loadVec(includes ,typeName ,constructors=None, methods=None):
        from dune.generator.generator import SimpleGenerator
    from dune.common.hashit import hashIt
    generator = SimpleGenerator("FieldVector", "Dune::Python")
    includes = includes + ["dune/python/common/fvector.hh"]
    typeHash = "fieldvector_" + hashIt(typeName)
    return generator.load(
        includes, typeName, typeHash,
        constructors, methods, bufferProtocol=True
    )


def FieldVector(values):
    """Construct a new FieldVector"""

    values = list(values)
    fv = "FieldVector_" + str(len(values))
    try:
        return globals()[fv](values)
    except KeyError:
        typeName = "Dune::FieldVector< double ," + str(len(values)) + " >"
        includes = []
        cls = _loadVec(includes, typeName).FieldVector
        setattr(cls, "_getitem", cls.__getitem__)
        setattr(cls, "__getitem__", _fieldVectorGetItem)
        globals().update({fv: cls})
    return globals()[fv](values)


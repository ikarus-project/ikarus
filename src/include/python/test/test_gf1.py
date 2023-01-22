# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import ikarus as iks
import numpy

if __name__ == "__main__":
    help(iks)
    assert iks.add(3,4)==7
    assert str(iks.scalarAffordances.noAffordance) == "scalarAffordances.noAffordance"
    E = numpy.array([[1, 7], [7, 4]])
    Evoigt= iks.to_voigt(E)
    assert Evoigt[0]==1
    assert Evoigt[1]==4
    assert Evoigt[2]==14

    E = numpy.array([[1, 7, 3],[7, 7, 3],[3, 3, 3]])
    Evoigt= iks.to_voigt(E)
    print(E)
    print(Evoigt)
    assert Evoigt[0]==1
    assert Evoigt[1]==7
    assert Evoigt[2]==3
    assert Evoigt[3]==6
    assert Evoigt[4]==6
    assert Evoigt[5]==14

    C = iks.planeStressLinearElasticMaterialTangent(1000,0.3)
    print(C)




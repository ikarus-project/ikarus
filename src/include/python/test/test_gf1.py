# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import ikarus
import numpy

if __name__ == "__main__":
    help(ikarus)
    assert ikarus.add(3,4)==7
    assert str(ikarus.scalarAffordances.noAffordance) == "scalarAffordances.noAffordance"
    E = numpy.array([[1, 7], [7, 4]])
    Evoigt= ikarus.to_voigt(E)
    assert Evoigt[0]==1
    assert Evoigt[1]==4
    assert Evoigt[2]==14




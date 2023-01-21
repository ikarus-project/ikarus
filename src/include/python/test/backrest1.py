# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import pickle, numpy
import dune.generator

def backup():
    from dune.grid import structuredGrid
    grid = structuredGrid([0,0],[1,1],[2,2])
    grid.hierarchicalGrid.globalRefine(2)
    a = numpy.array([1,2,3])

    pickle.dump([a,"hallo",grid.hierarchicalGrid,10], open("dumpA",'wb'))
    return grid

def restore():
    return pickle.load(open("dumpA","rb"))

class Test:
  def __init__(self,g,og):
      # note: only classes containing HGrids can be pickeled
      self.hg = g.hierarchicalGrid
      self.hog = og.hierarchicalGrid
  def run(self):
      return self.hg.leafView.size(0) == self.hog.leafView.size(0)



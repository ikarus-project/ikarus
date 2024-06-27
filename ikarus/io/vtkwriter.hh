// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vtkwriter.hh
 * \brief Ikarus VTK Writer for finite element results
 * \ingroup io
 *
 */

#include <dune/vtk/vtkwriter.hh>
#include <dune/vtk/writers/unstructuredgridwriter.hh>

template <typename GridView, bool structured = false, typename DC = void>
class VtkWroter
{
  // We are using the provided DataCollector, but if none was provided we are using the default ones from dune-vtk
  using DataCollector =
      std::conditional_t<std::is_same_v<DC, void>,
                         std::conditional_t<structured, typename Dune::Vtk::ContinuousDataCollector<GridView>,
                                            typename Dune::Vtk::YaspDataCollector<GridView>>,
                         DC>;

  // We are using a RectilinearGridWriter if structured is true
  using UnderlyingVTKWriter =
      std::conditional_t<structured, typename Dune::Vtk::RectilinearGridWriter<GridView, DataCollector>,
                         typename Dune::Vtk::UnstructuredGridWriter<GridView, DataCollector>>;

private:
  UnderlyingVTKWriter writer_;
};

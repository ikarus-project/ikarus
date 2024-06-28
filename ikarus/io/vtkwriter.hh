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

namespace Ikarus::Concepts {

// This has to be moved to simpleassemblers.hh
template <class AS>
concept IsAssembler = requires(AS as) {
  typename AS::GridView;
  as.finiteElements();
  as.requirement();
};

// taken from /dune/dune-vtk/dune/vtk/utility/concepts.hh
template <class DC>
concept IsDataCollector = requires(DC dc) {
  typename DC::GridView;
  dc.update();
  dc.numPoints();
  dc.numCells();
};
} // namespace Ikarus::Concepts

namespace Ikarus::Vtk {

template <typename AS, bool structured = false, typename DC = void>
requires(Concepts::IsAssembler<AS> and (std::is_same_v<DC, void> or Concepts::IsDataCollector<DC>))
class Writer
{
public:
  using Assembler     = AS;
  using GridView      = typename Assembler::GridView;
  using FERequirement = typename Assembler::FERequirement;
  using FEContainer   = typename Assembler::FEContainerType;
  using FEType        = typename std::remove_cvref_t<FEContainer>::value_type;

  // We are using the provided DataCollector, but if none was provided we are using the default ones from dune-vtk
  using DataCollector =
      std::conditional_t<std::is_same_v<DC, void>,
                         std::conditional_t<structured, typename Dune::Vtk::YaspDataCollector<GridView>,
                                            typename Dune::Vtk::ContinuousDataCollector<GridView>>,
                         DC>;

  // We are using a RectilinearGridWriter if structured is true
  using UnderlyingVTKWriter =
      std::conditional_t<structured, typename Dune::Vtk::RectilinearGridWriter<GridView, DataCollector>,
                         typename Dune::Vtk::UnstructuredGridWriter<GridView, DataCollector>>;

  template <class... Args>
  Writer(std::shared_ptr<Assembler> assembler, Args... args)
      : writer_(assembler->gridView(), args...),
        assembler_(assembler) {}

  template <class... Args>
  Writer(std::shared_ptr<Assembler> assembler, DataCollector& dc, Args... args)
      : writer_(dc, args...),
        assembler_(assembler) {}

  template <class... Args>
  Writer(std::shared_ptr<Assembler> assembler, DataCollector&& dc, Args... args)
      : writer_(std::move(dc), args...),
        assembler_(assembler) {}

  template <typename RF>
  void addResultFunction(RF&& resultFunction) {}

private:
  const FEContainer& finiteElements() { return assembler_->finiteElements(); }
  const FERequirement& requirement() { return assembler_->requirement(); }

  UnderlyingVTKWriter writer_;
  std::shared_ptr<Assembler> assembler_;
};

// Deduction guide

template <typename Assembler, class... Args>
requires(Concepts::IsAssembler<Assembler>)
Writer(std::shared_ptr<Assembler>, Args...) -> Writer<Assembler>;

template <typename DC, typename Assembler, class... Args, Dune::Vtk::IsDataCollector<DC> = true>
requires(Concepts::IsAssembler<Assembler>)
Writer(std::shared_ptr<Assembler>, DC&, Args...) -> Writer<Assembler, false, DC>;

template <typename DC, typename Assembler, class... Args, Dune::Vtk::IsDataCollector<DC> = true>
requires(Concepts::IsAssembler<Assembler>)
Writer(std::shared_ptr<Assembler>, DC&&, Args...) -> Writer<Assembler, false, DC>;

} // namespace Ikarus::Vtk

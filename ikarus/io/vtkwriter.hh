// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vtkwriter.hh
 * \brief Ikarus VTK Writer for finite element results
 * \ingroup io
 *
 */

#pragma once

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/vtk/vtkwriter.hh>
#include <dune/vtk/writers/unstructuredgridwriter.hh>

#include <ikarus/io/resultfunction.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus::Concepts {

// adapted from /dune/dune-vtk/dune/vtk/utility/concepts.hh
template <class DC>
concept IsDataCollector = requires(DC dc) {
  typename DC::GridView;
  dc.update();
  dc.numPoints();
  dc.numCells();
};
} // namespace Ikarus::Concepts

namespace Ikarus::Vtk {

namespace Impl {
  struct AsCellData
  {
  };
  struct AsPointData
  {
  };
  namespace Concepts {
    template <typename DT>
    concept DataType = std::is_same_v<DT, AsCellData> or std::is_same_v<DT, AsPointData>;
  }
} // namespace Impl

template <typename AS, bool structured = false, typename DC = void>
requires(Concepts::FlatAssembler<AS> and (std::is_same_v<DC, void> or Concepts::IsDataCollector<DC>))
class Writer
{
public:
  using Assembler     = AS;
  using GridView      = typename Assembler::GridView;
  using FERequirement = typename Assembler::FERequirement;
  using FEContainer   = typename Assembler::FEContainerType;
  using FEType        = typename std::remove_cvref_t<FEContainer>::value_type;

  // We are using the provided DataCollector, but if none was provided we are using the default ones from dune-vtk
  using DataCollector = std::conditional_t<std::is_same_v<DC, void>,
                                           std::conditional_t<structured, Dune::Vtk::YaspDataCollector<GridView>,
                                                              Dune::Vtk::ContinuousDataCollector<GridView>>,
                                           DC>;

  // We are using a RectilinearGridWriter if structured is true
  using UnderlyingVTKWriter = std::conditional_t<structured, Dune::Vtk::RectilinearGridWriter<GridView, DataCollector>,
                                                 Dune::Vtk::UnstructuredGridWriter<GridView, DataCollector>>;

  template <class... Args>
  explicit Writer(std::shared_ptr<Assembler> assembler, Args... args)
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

  // Sets the VTK file format
  void setFormat(Dune::Vtk::FormatTypes format) { writer_.setFormat(format); }

  /// Sets the global datatype used for coordinates and other global float values
  void setDatatype(Dune::Vtk::DataTypes datatype) { writer_.setDatatype(datatype); }

  /// Sets the integer type used in binary data headers
  void setHeadertype(Dune::Vtk::DataTypes datatype) { writer_.setHeadertype(datatype); }

  auto& writer() const { return writer_; }
  auto& writer() { return writer_; }

  template <typename RF, Impl::Concepts::DataType DT>
  void addResultFunction(RF&& resultFunction, DT /*type */) {
    if constexpr (std::is_same_v<DT, Impl::AsCellData>)
      writer_.addCellData(std::forward<RF>(resultFunction));
    else
      writer_.addPointData(std::forward<RF>(resultFunction));
  }

  template <template <typename, int, int> class RT, Impl::Concepts::DataType DT>
  requires(Concepts::ResultType<RT>)
  void addResult(DT type) {
    auto resFunction = makeResultVtkFunction<RT>(assembler_);
    addResultFunction(std::move(resFunction), type);
  }

  template <Impl::Concepts::DataType DT>
  void addAllResults(DT type) {
    using ResultTuple = typename FEType::SupportedResultTypes;

    Dune::Hybrid::forEach(ResultTuple(), [&]<typename RT>(RT i) { addResult<RT::template Rebind>(type); });
  }

  template <typename GF, Impl::Concepts::DataType DT>
  void addGridFunction(GF&& gridFunction, const std::string& name, size_t size, DT /* type */) {
    auto fieldInfo = Dune::Vtk::FieldInfo(name, size);
    if constexpr (std::is_same_v<DT, Impl::AsCellData>)
      writer_.addCellData(std::forward<GF>(gridFunction), fieldInfo);
    else
      writer_.addPointData(std::forward<GF>(gridFunction), fieldInfo);
  }

  template <typename GF, typename FieldInfo, Impl::Concepts::DataType DT>
  requires(std::same_as<FieldInfo, Dune::VTK::FieldInfo> or std::same_as<FieldInfo, Dune::Vtk::FieldInfo>)
  void addGridFunction(GF&& gridFunction, FieldInfo fieldInfo, DT /* type */) {
    if constexpr (std::is_same_v<DT, Impl::AsCellData>)
      writer_.addCellData(std::forward<GF>(gridFunction), fieldInfo);
    else
      writer_.addPointData(std::forward<GF>(gridFunction), fieldInfo);
  }

  template <int dim, typename R, typename Basis, Impl::Concepts::DataType DT>
  void addInterpolation(R&& vals, const Basis& basis, const std::string& name, DT type) {
    auto gridFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, dim>>(basis, std::forward<R>(vals));
    addGridFunction(std::move(gridFunction), name, dim, type);
  }

  // Other write functions are not explicitly exposed but can be access through writer(). ...
  std::string write(const std::string& fn, std::optional<std::string> dir = {}) const { return writer_.write(fn, dir); }

  // Legacy functions to add Cell and Point Data
  template <class... Args>
  void addCellData(Args&&... args) {
    writer_.addCellData(std::forward<Args>(args)...);
  }

  template <class... Args>
  void addPointData(Args&&... args) {
    writer_.addPointData(std::forward<Args>(args)...);
  }

private:
  const FEContainer& finiteElements() { return assembler_->finiteElements(); }
  const FERequirement& requirement() { return assembler_->requirement(); }

  UnderlyingVTKWriter writer_;
  std::shared_ptr<Assembler> assembler_;
};

// Helpers
inline auto asCellData() { return Impl::AsCellData{}; }

inline auto asPointData() { return Impl::AsPointData{}; }

// Deduction guide

template <typename Assembler, class... Args>
requires(Concepts::FlatAssembler<Assembler>)
Writer(std::shared_ptr<Assembler>, Args...) -> Writer<Assembler>;

template <typename DC, typename Assembler, class... Args, Dune::Vtk::IsDataCollector<DC> = true>
requires(Concepts::FlatAssembler<Assembler>)
Writer(std::shared_ptr<Assembler>, DC&, Args...) -> Writer<Assembler, false, DC>;

template <typename DC, typename Assembler, class... Args, Dune::Vtk::IsDataCollector<DC> = true>
requires(Concepts::FlatAssembler<Assembler>)
Writer(std::shared_ptr<Assembler>, DC&&, Args...) -> Writer<Assembler, false, DC>;

} // namespace Ikarus::Vtk

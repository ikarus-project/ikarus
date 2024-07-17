// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vtkwriter.hh
 * \brief Ikarus VTK Writer for finite element results
 * \ingroup io
 *
 */

#pragma once

#include <array>
#include <cstddef>
#include <tuple>

#include "dune/common/fvector.hh"
#include "dune/common/math.hh"
#include "dune/functions/functionspacebases/powerbasis.hh"
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/vtkwriter.hh>
#include <dune/vtk/writers/unstructuredgridwriter.hh>

#include <ikarus/io/resultfunction.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus::Concepts {} // namespace Ikarus::Concepts

// inline auto asCellData() { return Impl::AsCellData{}; }

// inline auto asPointData() { return Impl::AsPointData{}; }

namespace Ikarus::Vtk {
struct AsCellData
{
} asCellData;
struct AsPointData
{
} asPointData;

namespace Impl {

  template <typename GV, bool structured>
  struct DefaultDataCollector
  {
    using Type =
        std::conditional_t<structured, Dune::Vtk::YaspDataCollector<GV>, Dune::Vtk::ContinuousDataCollector<GV>>;
  };

  template <typename DC, bool structured>
  struct UnderlyingVTKWriter
  {
    using Type = std::conditional_t<structured, Dune::Vtk::RectilinearGridWriter<typename DC::GridView, DC>,
                                    Dune::Vtk::UnstructuredGridWriter<typename DC::GridView, DC>>;
  };

  namespace Concepts {
    template <typename DT>
    concept DataType = std::is_same_v<DT, AsCellData> or std::is_same_v<DT, AsPointData>;
  }

  template <typename Basis>
  constexpr auto dimBasis = []() {
    // Case 1 PowerBasis
    if constexpr (requires { Basis::PreBasis::children; })
      return Basis::PreBasis::children;
    // Case 2 SubSpaceBasis or Scalar Basis
    else
      return 1;
  }();

  template <typename Container>
  constexpr auto sizeOfContainer = []() {
    if constexpr (requires { Container::dimension; })
      return Container::dimension;
    else if constexpr (requires { std::tuple_size<Container>::value; })
      return std::tuple_size<Container>::value;
    else
      return 1ul;
  }();

} // namespace Impl

template <typename AS, bool structured = false,
          typename DC = Impl::DefaultDataCollector<typename AS::GridView, structured>::Type>
requires(Concepts::FlatAssembler<AS> and Concepts::DataCollector<DC>)
class Writer : public Impl::UnderlyingVTKWriter<DC, structured>::Type
{
public:
  using Assembler     = AS;
  using GridView      = typename Assembler::GridView;
  using FERequirement = typename Assembler::FERequirement;
  using FEContainer   = typename Assembler::FEContainerType;
  using FEType        = typename std::remove_cvref_t<FEContainer>::value_type;

  using DataCollector = DC;
  using Base          = Impl::UnderlyingVTKWriter<DC, structured>::Type;

  template <class... Args>
  explicit Writer(std::shared_ptr<Assembler> assembler, Args... args)
      : Base(assembler->gridView(), args...),
        assembler_(assembler) {}

  template <class... Args>
  Writer(std::shared_ptr<Assembler> assembler, DataCollector& dc, Args... args)
      : Base(dc, args...),
        assembler_(assembler) {}

  template <class... Args>
  Writer(std::shared_ptr<Assembler> assembler, DataCollector&& dc, Args... args)
      : Base(std::move(dc), args...),
        assembler_(assembler) {}

  template <typename RF, Impl::Concepts::DataType DT>
  void addResultFunction(RF&& resultFunction, DT type) {
    if constexpr (std::is_same_v<DT, AsCellData>)
      Base::addCellData(std::forward<RF>(resultFunction));
    else
      Base::addPointData(std::forward<RF>(resultFunction));
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

  template <typename Basis, typename Container = Dune::FieldVector<double, Impl::dimBasis<Basis>>,
            Impl::Concepts::DataType DT, typename R>
  void addInterpolation(R&& vals, const Basis& basis, const std::string& name, DT type) {
    auto gridFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Container>(basis, std::forward<R>(vals));
    auto fieldInfo    = Dune::Vtk::FieldInfo(name, Impl::sizeOfContainer<Container>);

    if constexpr (std::is_same_v<DT, AsCellData>)
      Base::addCellData(std::move(gridFunction), fieldInfo);
    else
      Base::addPointData(std::move(gridFunction), fieldInfo);
  }

private:
  std::shared_ptr<Assembler> assembler_;
};

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

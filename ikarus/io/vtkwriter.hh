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
#include <tuple>
#include <dune/common/fvector.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/vtkwriter.hh>
#include <dune/vtk/writers/unstructuredgridwriter.hh>

#include <ikarus/io/resultfunction.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus::Vtk {

struct AsCellData
{
} asCellData;
struct AsPointData
{
} asPointData;

namespace Impl {

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

template <typename AS, typename DC, typename Base>
requires(Concepts::FlatAssembler<AS> && Concepts::DataCollector<DC>)
struct Writer : public Base
{
public:
  using Assembler     = AS;
  using GridView      = typename Assembler::GridView;
  using FERequirement = typename Assembler::FERequirement;
  using FEContainer   = typename Assembler::FEContainer;
  using FEType        = typename std::remove_cvref_t<FEContainer>::value_type;

  using DataCollector = DC;
  using VTKWriter     = Base;

  template <class... Args>
  Writer(std::shared_ptr<AS> assembler, Args... args)
      : Base(assembler->gridView(), std::forward<Args>(args)...),
        assembler_(assembler) {}

  template <class... Args>
  Writer(std::shared_ptr<AS> assembler, DC& dc, Args... args)
      : Base(dc, std::forward<Args>(args)...),
        assembler_(assembler) {}

  template <class... Args>
  Writer(std::shared_ptr<AS> assembler, DC&& dc, Args... args)
      : Base(std::move(dc), std::forward<Args>(args)...),
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

/**
 * \brief Meta type to check whether a grid is structured, inherits from false_type
 *
 * \tparam G Grid type
 */
template <typename G>
struct IsStructured : std::false_type
{
};

/**
 * \brief Specialization of IsStructured for YASPGrids, inherits from true_type
 */
template <int dim, typename Coordinates>
struct IsStructured<Dune::YaspGrid<dim, Coordinates>> : std::true_type
{
};

/**
 * \brief Manages the default template parameter for the `Vtk::Writer`
 *
 * \tparam GV given GridView type
 */
template <typename GV>
requires Concepts::GridView<GV>
struct DefaultVTKWriterManager
{
  static constexpr bool isStructured = IsStructured<typename GV::Grid>::value;
  using DefaultDataCollector =
      std::conditional_t<isStructured, Dune::Vtk::YaspDataCollector<GV>, Dune::Vtk::ContinuousDataCollector<GV>>;

  template <typename DC = DefaultDataCollector>
  using DefaultVTKWriter = std::conditional_t<isStructured, Dune::Vtk::RectilinearGridWriter<typename DC::GridView, DC>,
                                              Dune::Vtk::UnstructuredGridWriter<typename DC::GridView, DC>>;
};

// Class template argument deduction guides for VTK::Writer

template <typename AS, class... Args>
requires(Ikarus::Concepts::FlatAssembler<AS>)
Writer(std::shared_ptr<AS>,
       Args...) -> Writer<AS, typename DefaultVTKWriterManager<typename AS::GridView>::DefaultDataCollector,
                          typename DefaultVTKWriterManager<typename AS::GridView>::template DefaultVTKWriter<>>;

template <typename AS, typename DC, class... Args, Dune::Vtk::IsDataCollector<DC> = true>
requires(Ikarus::Concepts::FlatAssembler<AS>)
Writer(std::shared_ptr<AS>, DC&, Args...)
    -> Writer<AS, DC, typename DefaultVTKWriterManager<typename AS::GridView>::template DefaultVTKWriter<DC>>;

template <typename AS, typename DC, class... Args, Dune::Vtk::IsDataCollector<DC> = true>
requires(Ikarus::Concepts::FlatAssembler<AS>)
Writer(std::shared_ptr<AS>, DC&&, Args...)
    -> Writer<AS, DC, typename DefaultVTKWriterManager<typename AS::GridView>::template DefaultVTKWriter<DC>>;

} // namespace Ikarus::Vtk
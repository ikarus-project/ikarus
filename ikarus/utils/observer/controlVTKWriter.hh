// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "observer.hh"
#include "observerMessages.hh"

#include <string>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <spdlog/spdlog.h>
#include <ikarus/io/shell2DDataCollector.hh>

template <typename Basis, typename SolutionVectorType>  // Check basis
class ControlSubsamplingVertexVTKWriter : public IObserver<ControlMessages> {

  using DC = Dune::Vtk::Shell2DDataCollector<typename Basis::GridView>;
public:
  template< typename FunctionType>
  ControlSubsamplingVertexVTKWriter(const Basis& p_basis, const SolutionVectorType& sol,FunctionType&& p_func, int refinementLevels = 0)
      : basis{&p_basis}, vtkWriter(std::make_shared<DC>(p_basis.gridView(),Dune::refinementLevels(refinementLevels)),Dune::Vtk::FormatTypes::ASCII,Dune::Vtk::DataTypes::FLOAT64), solution{&sol}, func{p_func} {}


  auto setFileNamePrefix(std::string&& p_name) { prefixString = std::move(p_name); }

  void updateImpl(ControlMessages message) override {
    switch (message) {
      case ControlMessages::SOLUTION_CHANGED: {
        func(vtkWriter, *basis, *solution, prefixString,step++);

      } break;
      default:
        break;  //   default: do nothing when notified
    }
  }

  using IObserver::updateImpl;
  void updateImpl(ControlMessages, double) override {}
  void updateImpl(ControlMessages, const Eigen::VectorXd&) override {}

private:
  Basis const* basis;
  using VtkWriter = Dune::VtkUnstructuredGridWriter<typename Basis::GridView,DC>;
  VtkWriter vtkWriter;
  SolutionVectorType const* solution;
  int step{0};
  std::function<void(VtkWriter&,const Basis&, const SolutionVectorType&, std::string&,int)> func;
  std::string prefixString{};
};

//
// Created by lex on 14/12/2021.
//

#pragma once
#include <string>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "spdlog/spdlog.h"

#include "ikarus/utils/Observer/observer.h"
#include "ikarus/utils/Observer/observerMessages.h"

template <typename Basis>  // Check basis
class ControlSubsamplingVertexVTKWriter : public IObserver<ControlMessages> {
public:
  ControlSubsamplingVertexVTKWriter(const Basis& p_basis, const Eigen::VectorXd& sol, int refinementLevels)
      : basis{&p_basis}, vtkWriter(p_basis.gridView(), Dune::refinementLevels(refinementLevels)), solution{&sol} {}

  void setVertexSolutionName(std::string&& p_name) { solutionName = std::move(p_name); }

  void setFileNamePrefix(std::string&& p_name) { prefixString = std::move(p_name); }

  void updateImpl(ControlMessages message) override {
    switch (message) {
      case ControlMessages::SOLUTION_CHANGED: {
        auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<
            Dune::FieldVector<double, Basis::LocalView::Tree::CHILDREN>>(*basis, *solution);
        vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo(solutionName, Dune::VTK::FieldInfo::Type::vector,
                                                           Basis::LocalView::Tree::CHILDREN));
        vtkWriter.write(prefixString + std::to_string(step++));
      } break;
      default:
        break;  //   default: do nothing when notified
    }
  }

  void updateImpl(ControlMessages, double) override {}
  void updateImpl(ControlMessages, const Eigen::VectorXd&) override {}

private:
  Basis const* basis;
  Dune::SubsamplingVTKWriter<typename Basis::GridView> vtkWriter;
  Eigen::VectorXd const* solution;
  int step{0};
  std::string solutionName{"displacement"};
  std::string prefixString{};
};

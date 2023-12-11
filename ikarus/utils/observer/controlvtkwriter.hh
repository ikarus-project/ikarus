// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "observer.hh"
#include "observermessages.hh"

#include <string>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <spdlog/spdlog.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {

  template <typename Basis>  // Check basis
  class ControlSubsamplingVertexVTKWriter : public IObserver<ControlMessages> {
    static constexpr int components = Basis::LocalView::Tree::degree() == 0 ? 1 : Basis::LocalView::Tree::degree();

  public:
    ControlSubsamplingVertexVTKWriter(const Basis& p_basis, const Eigen::VectorXd& sol, int refinementLevels = 0)
        : basis{&p_basis}, vtkWriter(p_basis.gridView(), Dune::refinementLevels(refinementLevels)), solution{&sol} {}

    auto setFieldInfo(std::string&& name, Dune::VTK::FieldInfo::Type type, std::size_t size,
                      Dune::VTK::Precision prec = Dune::VTK::Precision::float32) {
      fieldInfo      = Dune::VTK::FieldInfo(std::move(name), type, size, prec);
      isFieldInfoSet = true;
    }

    auto setFileNamePrefix(std::string&& p_name) { prefixString = std::move(p_name); }

    void updateImpl(ControlMessages message) override {
      assert(isFieldInfoSet && "You need to call setFieldInfo first!");
      switch (message) {
        case ControlMessages::SOLUTION_CHANGED: {
          auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, components>>(
              *basis, *solution);
          vtkWriter.addVertexData(disp, fieldInfo);
          vtkWriter.write(prefixString + std::to_string(step++));
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
    Dune::SubsamplingVTKWriter<typename Basis::GridView> vtkWriter;
    Eigen::VectorXd const* solution;
    int step{0};
    Dune::VTK::FieldInfo fieldInfo{"Default", Dune::VTK::FieldInfo::Type::scalar, 1};
    std::string prefixString{};
    bool isFieldInfoSet{false};
  };
}  // namespace Ikarus
#pragma GCC diagnostic pop

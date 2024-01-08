// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlvtkwriter.hh
 * \brief Observer implementation for writing vtk files when notified
 */

#pragma once
#include "observer.hh"
#include "observermessages.hh"

#include <string>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {

  /**
   * @brief ControlSubsamplingVertexVTKWriter class for writing VTK files with subsampling based on control messages.
   *
   * \details It inherits from the IObserver class and is specifically designed for handling SOLUTION_CHANGED messages.
   *
   * @tparam Basis The type of the grid basis.
   */
  template <typename Basis>
  class ControlSubsamplingVertexVTKWriter : public IObserver<ControlMessages> {
    static constexpr int components = Basis::LocalView::Tree::degree() == 0 ? 1 : Basis::LocalView::Tree::degree();

  public:
    /**
     * @brief Constructor for ControlSubsamplingVertexVTKWriter.
     *
     * Initializes the VTK writer with the provided basis, solution, and refinement levels.
     *
     * @param p_basis The grid basis.
     * @param sol The solution vector.
     * @param refinementLevels The refinement levels for subsampling.
     */
    ControlSubsamplingVertexVTKWriter(const Basis& p_basis, const Eigen::VectorXd& sol, int refinementLevels = 0)
        : basis{&p_basis}, vtkWriter(p_basis.gridView(), Dune::refinementLevels(refinementLevels)), solution{&sol} {}

    /**
     * @brief Set field information for the VTK file.
     *
     * @param name The name of the field.
     * @param type The type of the field.
     * @param size The size of the field.
     * @param prec The precision of the field.
     * @return The field information.
     */
    auto setFieldInfo(std::string&& name, Dune::VTK::FieldInfo::Type type, std::size_t size,
                      Dune::VTK::Precision prec = Dune::VTK::Precision::float32) {
      fieldInfo      = Dune::VTK::FieldInfo(std::move(name), type, size, prec);
      isFieldInfoSet = true;
    }

    /**
     * @brief Set the file name prefix for VTK files.
     *
     * @param p_name The file name prefix.
     */
    auto setFileNamePrefix(std::string&& p_name) { prefixString = std::move(p_name); }

    /**
     * @brief Implementation of the update method.
     *
     * This method is called upon receiving a SOLUTION_CHANGED control message.
     * It writes VTK files with subsampling based on the provided field information.
     *
     * @param message The received control message.
     */
    void updateImpl(ControlMessages message) final {
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

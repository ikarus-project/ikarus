// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file dynamics.hh
 * \brief Helper for
 */

#pragma once
#include <matplot/matplot.h>

#include <dune/vtk/pvdwriter.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/solver/eigenvaluesolver/generaleigensolver.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

template <typename FEC, typename DV>
struct ModalAnalysis
{
  using Assembler     = SparseFlatAssembler<FEC, DV>;
  using FERequirement = typename Assembler::FERequirement;
  using LumpedAssembler =
      AssemblerManipulator<Assembler, Ikarus::Impl::AssemblerInterfaceHelper<ScalarAssembler, ScalarManipulator>,
                           Ikarus::Impl::AssemblerInterfaceHelper<VectorAssembler, VectorManipulator>,
                           Ikarus::Impl::AssemblerInterfaceHelper<MatrixAssembler, MatrixManipulator>>;
  using Solver = GeneralSymEigenSolver<EigenSolverTypeTag::Spectra, Eigen::SparseMatrix<double>>;

  template <typename FES>
  ModalAnalysis(FES&& fes, const DV& dv)
      : stiffAssembler_(makeSparseFlatAssembler(std::forward<FES>(fes), dv)),
        massAssembler_(makeSparseFlatAssembler(std::forward<FES>(fes), dv)) {
    d_.setZero(dv.basis().size());
    req_.insertGlobalSolution(d_);

    stiffAssembler_->bind(req_, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);
    massAssembler_->bind(req_, Ikarus::AffordanceCollections::dynamics, Ikarus::DBCOption::Reduced);
    lumpedMassAssembler_ = makeAssemblerManipulator(*massAssembler_);
  }

  template <typename LumpingScheme>
  void registerLumpingScheme(LumpingScheme ls = LumpingScheme{}) {
    lumpedMassAssembler_->unbindAllMatrixFunctions();
    lumpedMassAssembler_->bind(ls);
  }

  bool compute() {
    solver_.emplace(stiffAssembler_, lumpedMassAssembler_);
    return solver_->compute();
  }

  Eigen::VectorXd angularFrequencies() { return eigenvalues().cwiseSqrt().eval(); }

  Eigen::VectorXd naturalFrequencies() { return angularFrequencies() / (2 * std::numbers::pi); }

  const Eigen::VectorXd& eigenvalues() const { return solver_->eigenvalues(); }

  void plotModalSpectrum() {
    using namespace matplot;
    auto freq = naturalFrequencies();
    std::vector<double> frequencies(freq.data(), freq.data() + freq.size());

    auto modeNumbersView =
        std::ranges::iota_view{1ul, frequencies.size() + 1} |
        std::ranges::views::transform([&](size_t i) { return static_cast<double>(i) / frequencies.size(); });
    std::vector modeNumbers(modeNumbersView.begin(), modeNumbersView.end());

    auto fig = figure(true);
    auto ax  = fig->add_axes();
    ax->plot(modeNumbers, frequencies)->line_width(1).color("b");
    ax->xlabel("Mode Number");
    ax->ylabel("Frequency (Hz)");
    ax->title("Modal Analysis Spectrum");
    ax->grid(true);

    // Show the figure
    matplot::show();
  }

private:
  std::shared_ptr<Assembler> stiffAssembler_;
  std::shared_ptr<Assembler> massAssembler_;
  std::shared_ptr<LumpedAssembler> lumpedMassAssembler_;
  FERequirement req_{};
  typename FERequirement::SolutionVectorType d_;
  std::optional<Solver> solver_{};
};

template <typename FEC, typename DV>
ModalAnalysis(FEC&&, const DV&) -> ModalAnalysis<FEC, DV>;

} // namespace Ikarus::Dynamics

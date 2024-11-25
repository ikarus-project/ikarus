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
#include <ikarus/utils/dynamics/dynamicshelpers.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

MAKE_ENUM(ModalAnalysisResultType, squaredAngularFrequency, angularFrequency, naturalFrequency);

template <typename FEC, typename DV>
struct ModalAnalysis
{
  using Assembler     = SparseFlatAssembler<FEC, DV>;
  using MatrixType    = typename Assembler::MatrixType;
  using FERequirement = typename Assembler::FERequirement;
  using FEContainer   = FEC;

  // static_assert(MatrixType::IsRowMajor, "Rowmajor");

  using LumpedAssembler =
      AssemblerManipulator<Assembler, Ikarus::Impl::AssemblerInterfaceHelper<ScalarAssembler, ScalarManipulator>,
                           Ikarus::Impl::AssemblerInterfaceHelper<VectorAssembler, VectorManipulator>,
                           Ikarus::Impl::AssemblerInterfaceHelper<MatrixAssembler, MatrixManipulator>>;
  using Solver = GeneralSymEigenSolver<EigenSolverTypeTag::Spectra, MatrixType>;

  template <typename FES>
  ModalAnalysis(FES&& fes, const DV& dv)
      : fes_(std::forward<FES>(fes)),
        d_(dv.basis().size()),
        dRef_(d_),
        stiffAssembler_(std::make_shared<Assembler>(fes_, dv)),
        massAssembler_(std::make_shared<Assembler>(fes_, dv)) {
    if constexpr (std::remove_cvref_t<decltype(fes.front())>::Traits::useEigenRef)
      req_.insertGlobalSolution(dRef_);
    else
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

  Eigen::VectorXd angularFrequencies() {
    assertCompute();
    return squaredAngularFrequencies().cwiseSqrt().eval();
  }

  Eigen::VectorXd naturalFrequencies() {
    assertCompute();
    return angularFrequencies() / (2 * std::numbers::pi);
  }

  const Eigen::VectorXd& squaredAngularFrequencies() const {
    assertCompute();
    return solver_->eigenvalues();
  }

  auto frequencies(ModalAnalysisResultType rt) {
    if (rt == ModalAnalysisResultType::angularFrequency)
      return angularFrequencies();
    if (rt == ModalAnalysisResultType::naturalFrequency)
      return naturalFrequencies();
    if (rt == ModalAnalysisResultType::squaredAngularFrequency)
      return squaredAngularFrequencies();
    DUNE_THROW(Dune::NotImplemented, "Requested result not implemented");
  }

  const Eigen::MatrixXd& eigenmodes() const { return solver_->eigenvectors(); }

  void plotModalSpectrum(ModalAnalysisResultType resultType = ModalAnalysisResultType::angularFrequency,
                         bool normalizeModeNumber           = false) {
    assertCompute();
    using namespace matplot;
    auto freq = frequencies(resultType);

    auto modeNumbersView = std::ranges::iota_view{1l, nev() + 1} | std::ranges::views::transform([&](auto i) {
                             if (normalizeModeNumber)
                               return static_cast<double>(i) / nev();
                             return static_cast<double>(i);
                           });

    std::vector<double> modeNumbers(modeNumbersView.begin(), modeNumbersView.end());

    auto fig = figure(true);
    auto ax  = fig->add_axes();
    ax->plot(modeNumbers, freq)->line_width(1).color("b");
    ax->xlabel("Mode Number");
    ax->ylabel("Frequency (Hz)");
    ax->title("Modal Analysis Spectrum");
    ax->grid(true);

    // Show the figure
    matplot::show();
  }

  void writeEigenModes(const std::string& filename, std::optional<Eigen::Index> nev_ = std::nullopt) const {
    assertCompute();
    writeEigenmodesToPVD(solver_.value(), stiffAssembler_, filename, nev_);
  }

  auto nev() const { return solver_->nev(); }

private:
  FEContainer fes_;
  FERequirement req_{};
  Eigen::VectorXd d_;
  Eigen::Ref<Eigen::VectorXd> dRef_;
  std::shared_ptr<Assembler> stiffAssembler_;
  std::shared_ptr<Assembler> massAssembler_;
  std::shared_ptr<LumpedAssembler> lumpedMassAssembler_;

  std::optional<Solver> solver_{};

  void assertCompute() const {
    if (not solver_)
      DUNE_THROW(Dune::IOError, "Eigenvalues and -vectors not yet computed, please call compute() first");
  }
};

template <typename FEC, typename DV>
ModalAnalysis(FEC&&, const DV&) -> ModalAnalysis<FEC, DV>;

} // namespace Ikarus::Dynamics

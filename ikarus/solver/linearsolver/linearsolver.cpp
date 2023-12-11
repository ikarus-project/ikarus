// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "linearsolver.hh"

#include <Eigen/CholmodSupport>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/UmfPackSupport>
// #include <Eigen/SuperLUSupport>

namespace Ikarus {

  template <typename ScalarType>
  LinearSolverTemplate<ScalarType>::LinearSolverTemplate(const SolverTypeTag& p_solverTypeTag)
      : solverTypeTag{p_solverTypeTag} {
    using namespace Eigen;

    switch (solverTypeTag) {
      case SolverTypeTag::si_ConjugateGradient:
        solverimpl = std::make_unique<SolverImpl<ConjugateGradient<SparseMatrixType, Lower | Upper>>>();
        break;
      case SolverTypeTag::si_LeastSquaresConjugateGradient:
        solverimpl = std::make_unique<SolverImpl<LeastSquaresConjugateGradient<SparseMatrixType>>>();
        break;
      case SolverTypeTag::si_BiCGSTAB:
        solverimpl = std::make_unique<SolverImpl<BiCGSTAB<SparseMatrixType>>>();
        break;
      case SolverTypeTag::sd_SimplicialLLT:
        solverimpl = std::make_unique<SolverImpl<SimplicialLLT<SparseMatrixType>>>();
        break;
      case SolverTypeTag::sd_SimplicialLDLT:
        solverimpl = std::make_unique<SolverImpl<SimplicialLDLT<SparseMatrixType>>>();
        break;
      case SolverTypeTag::sd_SparseLU:
        solverimpl = std::make_unique<SolverImpl<SparseLU<SparseMatrixType>>>();
        break;
      case SolverTypeTag::sd_SparseQR:
        solverimpl = std::make_unique<SolverImpl<SparseQR<SparseMatrixType, COLAMDOrdering<int>>>>();
        break;
      case SolverTypeTag::sd_CholmodSupernodalLLT:
        solverimpl = std::make_unique<SolverImpl<CholmodSupernodalLLT<SparseMatrixType>>>();
        break;
      case SolverTypeTag::sd_UmfPackLU:
        solverimpl = std::make_unique<SolverImpl<UmfPackLU<SparseMatrixType>>>();
        break;
      case SolverTypeTag::sd_SuperLU:
        DUNE_THROW(Dune::NotImplemented, "Not implemented yet.");
        break;
        // Dense Solver
      case SolverTypeTag::d_PartialPivLU:
        solverimpl = std::make_unique<SolverImpl<PartialPivLU<DenseMatrixType>>>();
        break;
      case SolverTypeTag::d_FullPivLU:
        solverimpl = std::make_unique<SolverImpl<FullPivLU<DenseMatrixType>>>();
        break;
      case SolverTypeTag::d_HouseholderQR:
        solverimpl = std::make_unique<SolverImpl<HouseholderQR<DenseMatrixType>>>();
        break;
      case SolverTypeTag::d_ColPivHouseholderQR:
        solverimpl = std::make_unique<SolverImpl<ColPivHouseholderQR<DenseMatrixType>>>();
        break;
      case SolverTypeTag::d_FullPivHouseholderQR:
        solverimpl = std::make_unique<SolverImpl<FullPivHouseholderQR<DenseMatrixType>>>();
        break;
      case SolverTypeTag::d_CompleteOrthogonalDecomposition:
        solverimpl = std::make_unique<SolverImpl<CompleteOrthogonalDecomposition<DenseMatrixType>>>();
        break;
      case SolverTypeTag::d_LLT:
        solverimpl = std::make_unique<SolverImpl<LLT<DenseMatrixType>>>();
        break;
      case SolverTypeTag::d_LDLT:
        solverimpl = std::make_unique<SolverImpl<LDLT<DenseMatrixType>>>();
        break;
      case SolverTypeTag::none:
        break;
      default:
        DUNE_THROW(Dune::NotImplemented, "Your requested solver does not work with this interface class");
    }
  }
  template class LinearSolverTemplate<double>;

}  // namespace Ikarus

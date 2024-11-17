// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file generaleigensolver.hh
 * \brief Helper for dune-functions
 */

#pragma once

#include <optional>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>

#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

MAKE_ENUM(EigenSolverTypeTag, Spectra, Eigen);
MAKE_ENUM(MatrixTypeTag, Dense, Sparse);

template <EigenSolverTypeTag SolverType, MatrixTypeTag matrixType, typename ScalarType>
struct GeneralSymEigenSolver
{
};
/**
 * \brief 
 * 
 * \tparam matrixType 
 * \tparam ST
 * \ingroup Dynamics 
 */
template <MatrixTypeTag matrixType, typename ST>
struct GeneralSymEigenSolver<EigenSolverTypeTag::Spectra, matrixType, ST>
{
  using ScalarType              = ST;
  static constexpr bool isDense = matrixType == MatrixTypeTag::Dense;

  using MatrixType = std::conditional_t<isDense, Eigen::MatrixX<ScalarType>, Eigen::SparseMatrix<ScalarType>>;
  using ProductType =
      std::conditional_t<isDense, Spectra::DenseSymMatProd<ScalarType>, Spectra::SparseSymMatProd<ScalarType>>;
  using CholeskyType =
      std::conditional_t<isDense, Spectra::DenseCholesky<ScalarType>, Spectra::SparseCholesky<ScalarType>>;

  using SolverType = Spectra::SymGEigsSolver<ProductType, CholeskyType, Spectra::GEigsMode::Cholesky>;

  /**
   * \brief
   *
   * \tparam MATA
   * \tparam MATB
   */
  template <typename MATA, typename MATB>
  requires(std::convertible_to<MATA, MatrixType>, std::convertible_to<MATB, MatrixType>)
  GeneralSymEigenSolver(MATA&& A, MATB&& B, Eigen::Index nev)
      : nev_(nev),
        aOP_(std::forward<MATA>(A)),
        bOP_(std::forward<MATB>(B)),
        solver_(aOP_, bOP_, nev, 2 * nev <= A.cols() ? 2 * nev : A.cols()) {
    if (A.cols() != B.cols())
      DUNE_THROW(Dune::IOError, "GeneralSymEigenSolver: The passed matrices should have the same size");
  }

  template <Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
  GeneralSymEigenSolver(const std::shared_ptr<AssemblerA> assemblerA, const std::shared_ptr<AssemblerB>& assemblerB,
                        Eigen::Index nev)
      : GeneralSymEigenSolver(assemblerA->matrix(), assemblerB->matrix(), nev) {
    if (not(assemblerA->dBCOption() == DBCOption::Reduced && assemblerB->dBCOption() == DBCOption::Reduced))
      DUNE_THROW(Dune::IOError, "GeneralSymEigenSolver: The passed assembler should both have DBCOption::Reduced");
  }

  /**
   * \brief Starts the computation of the eigenvalue solver
   *
   * \param tolerance given tolerance for iterative eigenvalue solving (default: 1e-10)
   * \param maxit givenn maximum iterations for eigenvalue solving (default 1000)
   * \return true solving was successfull
   * \return false solving was not successfull
   */
  bool compute(ScalarType tolerance = 1e-10, Eigen::Index maxit = 1000) {
    solver_.init();
    solver_.compute(Spectra::SortRule::SmallestAlge, 1000, 1e-10, Spectra::SortRule::SmallestAlge);

    computed_ = solver_.info() == Spectra::CompInfo::Successful;
    return computed_;
  }

  /**
   * \brief Returns the eigenvalues of the gerneral eigenvalue problem
   *
   * \param angularFrequency if true, obtained eigenvalues are returned square-rooted
   * \return Eigen::VectorXd vector of eigenvalues
   */
  Eigen::VectorXd eigenvalues(bool angularFrequency = false) {
    computeIfNeeded();
    if (angularFrequency)
      return solver_.eigenvalues().array().sqrt().eval();
    return solver_.eigenvalues();
  }

  /**
   * \brief Returns the eigenvectors of the gerneral eigenvalue problem
   *
   * \param _nev optionally specify how many eigenvectors are requested
   * \return auto matrix with the eigevectors as columns
   */
  auto eigenvectors(std::optional<Eigen::Index> _nev = std::nullopt) {
    computeIfNeeded();
    return solver_.eigenvectors(_nev.value_or(nev_));
  }

private:
  Eigen::Index nev_;
  ProductType aOP_;
  CholeskyType bOP_;
  SolverType solver_;
  bool computed_{};

  void computeIfNeeded() {
    if (not computed_)
      compute();
  }
};

template <MatrixTypeTag matrixType, typename ST>
struct GeneralSymEigenSolver<EigenSolverTypeTag::Eigen, matrixType, ST>
{
};

} // namespace Ikarus::Dynamics

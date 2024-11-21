// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file generaleigensolver.hh
 * \brief Helper for dune-functions
 */

#pragma once

#include <optional>
#include <utility>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>

#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus {

MAKE_ENUM(EigenSolverTypeTag, Spectra, Eigen);

template <EigenSolverTypeTag SolverType, Concepts::DenseOrSparseEigenMatrix MT>
struct GeneralSymEigenSolver
{
};

template <Concepts::DenseOrSparseEigenMatrix MT>
struct GeneralSymEigenSolver<EigenSolverTypeTag::Spectra, MT>
{
  using ScalarType              = typename MT::Scalar;
  static constexpr bool isDense = Concepts::EigenMatrix<MT>;

  using MatrixType = MT;
  using ProductType =
      std::conditional_t<isDense, Spectra::DenseSymMatProd<ScalarType>, Spectra::SparseSymMatProd<ScalarType>>;
  using CholeskyType =
      std::conditional_t<isDense, Spectra::DenseCholesky<ScalarType>, Spectra::SparseCholesky<ScalarType>>;

  using SolverType = Spectra::SymGEigsSolver<ProductType, CholeskyType, Spectra::GEigsMode::Cholesky>;

  template <typename MATA, typename MATB>
  requires(Concepts::DenseOrSparseEigenMatrix<std::remove_cvref_t<MATA>>)
  GeneralSymEigenSolver(MATA&& A, MATB&& B)
      : nev_(A.rows()),
        nevsPartition_(static_cast<Eigen::Index>(std::ceil(nev_ / 2)),
                       static_cast<Eigen::Index>(nev_ - std::ceil(nev_ / 2))),
        aOP_(std::forward<MATA>(A)),
        bOP_(std::forward<MATB>(B)),
        solverSmallest_(aOP_, bOP_, nevsPartition_.first, nevsPartition_.second * 2),
        solverGreatest_(aOP_, bOP_, nevsPartition_.second, nevsPartition_.second * 2) {
    if (A.cols() != B.cols())
      DUNE_THROW(Dune::IOError, "GeneralSymEigenSolver: The passed matrices should have the same size");
    eigenvalues_.resize(nev_);
    eigenvectors_.resize(A.rows(), nev_);
  }

  template <Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
  GeneralSymEigenSolver(const std::shared_ptr<AssemblerA>& assemblerA, const std::shared_ptr<AssemblerB>& assemblerB)
      : GeneralSymEigenSolver(assemblerA->matrix(), assemblerB->matrix()) {
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
    solverSmallest_.init();
    solverSmallest_.compute(Spectra::SortRule::SmallestAlge, 1000, 1e-10, Spectra::SortRule::SmallestAlge);

    solverGreatest_.init();
    solverGreatest_.compute(Spectra::SortRule::LargestAlge, 1000, 1e-10, Spectra::SortRule::SmallestAlge);

    eigenvalues_.head(nevsPartition_.first)  = solverSmallest_.eigenvalues();
    eigenvalues_.tail(nevsPartition_.second) = solverGreatest_.eigenvalues();

    eigenvectors_.leftCols(nevsPartition_.first)   = solverSmallest_.eigenvectors();
    eigenvectors_.rightCols(nevsPartition_.second) = solverGreatest_.eigenvectors();

    computed_ = solverSmallest_.info() == Spectra::CompInfo::Successful and
                solverGreatest_.info() == Spectra::CompInfo::Successful;

    return computed_;
  }

  /**
   * \brief Returns the eigenvalues of the gerneral eigenvalue problem
   *
   * \return Eigen::VectorXd vector of eigenvalues
   */
  auto& eigenvalues() const {
    assertCompute();
    return eigenvalues_;
  }

  /**
   * \brief Returns the eigenvectors of the gerneral eigenvalue problem
   *
   * \param _nev optionally specify how many eigenvectors are requested
   * \return auto matrix with the eigevectors as columns
   */
  auto& eigenvectors() const {
    assertCompute();
    return eigenvectors_;
  }

  Eigen::Index nev() const { return nev_; }

private:
  Eigen::Index nev_;
  std::pair<Eigen::Index, Eigen::Index> nevsPartition_;
  ProductType aOP_;
  CholeskyType bOP_;
  SolverType solverSmallest_;
  SolverType solverGreatest_;

  Eigen::VectorXd eigenvalues_;
  Eigen::MatrixXd eigenvectors_;

  bool computed_{};

  void assertCompute() const {
    if (not computed_)
      DUNE_THROW(Dune::IOError, "Eigenvalues and -vectors not yet computed, please call compute() first");
  }
};

template <Concepts::EigenMatrix MT>
struct GeneralSymEigenSolver<EigenSolverTypeTag::Eigen, MT>
{
  using ScalarType = typename MT::Scalar;
  using MatrixType = MT;
  using SolverType = Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType>;

  template <typename MATA, typename MATB>
  requires(Concepts::EigenMatrix<std::remove_cvref_t<MATA>>)
  GeneralSymEigenSolver(MATA&& A, MATB&& B)
      : matA_(A),
        matB_(B),
        solver_(A.size()) {
    if (A.cols() != B.cols())
      DUNE_THROW(Dune::IOError, "GeneralSymEigenSolver: The passed matrices should have the same size");
  }

  template <Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
  GeneralSymEigenSolver(const std::shared_ptr<AssemblerA> assemblerA, const std::shared_ptr<AssemblerB>& assemblerB)
      : GeneralSymEigenSolver(assemblerA->matrix(), assemblerB->matrix()) {
    if (not(assemblerA->dBCOption() == DBCOption::Reduced && assemblerB->dBCOption() == DBCOption::Reduced))
      DUNE_THROW(Dune::IOError, "GeneralSymEigenSolver: The passed assembler should both have DBCOption::Reduced");
  }

  /**
   * \brief Starts the computation of the eigenvalue solver
   *
   * \param options defaults to Eigen::ComputeEigenvectors, can be set to Eigen::EigenvaluesOnly, accessing eigenvectors
   * in that case will results in an error
   * \return true solving was successfull
   * \return false solving was not successfull
   */
  bool compute(int options = Eigen::ComputeEigenvectors) {
    solver_.compute(matA_, matB_, options);
    computed_ = true;
    return true; // SelfAdjointEigenSolver will always be successfull if prerequisits are met
  }

  /**
   * \brief Returns the eigenvalues of the gerneral eigenvalue problem
   *
   * \return Reference to the vector of eigenvalues
   */
  auto& eigenvalues() const {
    assertCompute();
    return solver_.eigenvalues();
  }

  /**
   * \brief Returns the eigenvectors of the gerneral eigenvalue problem
   *
   * \param _nev optionally specify how many eigenvectors are requested
   * \return Reference to the matrix with the eigevectors as columns
   */
  auto& eigenvectors() const {
    assertCompute();
    return solver_.eigenvectors();
  }

private:
  MatrixType matA_;
  MatrixType matB_;
  SolverType solver_;
  bool computed_{};

  void assertCompute() const {
    if (not computed_)
      DUNE_THROW(Dune::IOError, "Eigenvalues and -vectors not yet computed, please call compute() first");
  }
};

template <EigenSolverTypeTag tag, Concepts::FlatAssembler AS1, Concepts::FlatAssembler AS2>
requires(std::same_as<typename AS1::MatrixType, typename AS2::MatrixType> &&
         not(tag == EigenSolverTypeTag::Eigen && Concepts::SparseEigenMatrix<typename AS1::MatrixType>))
auto makeGeneralSymEigenSolver(const std::shared_ptr<AS1>& as1, const std::shared_ptr<AS2> as2) {
  using MatrixType        = typename AS1::MatrixType;
  constexpr auto isSparse = Concepts::SparseEigenMatrix<MatrixType>;
  using SolverType =
      std::conditional_t<isSparse, GeneralSymEigenSolver<tag, MatrixType>, GeneralSymEigenSolver<tag, MatrixType>>;

  return SolverType{as1, as2};
}

template <Concepts::DenseOrSparseEigenMatrix MT>
struct PartialGeneralSymEigenSolver
{
  using ScalarType              = typename MT::Scalar;
  static constexpr bool isDense = Concepts::EigenMatrix<MT>;

  using MatrixType = MT;
  using ProductType =
      std::conditional_t<isDense, Spectra::DenseSymMatProd<ScalarType>, Spectra::SparseSymMatProd<ScalarType>>;
  using CholeskyType =
      std::conditional_t<isDense, Spectra::DenseCholesky<ScalarType>, Spectra::SparseCholesky<ScalarType>>;

  using SolverType = Spectra::SymGEigsSolver<ProductType, CholeskyType, Spectra::GEigsMode::Cholesky>;

  template <typename MATA, typename MATB>
  requires(Concepts::DenseOrSparseEigenMatrix<std::remove_cvref_t<MATA>>)
  PartialGeneralSymEigenSolver(MATA&& A, MATB&& B, Eigen::Index nev)
      : nev_(nev),
        aOP_(std::forward<MATA>(A)),
        bOP_(std::forward<MATB>(B)),
        solver_(aOP_, bOP_, nev, 2 * nev <= A.cols() ? 2 * nev : A.cols()) {
    if (A.cols() != B.cols())
      DUNE_THROW(Dune::IOError, "GeneralSymEigenSolver: The passed matrices should have the same size");
  }

  template <Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
  PartialGeneralSymEigenSolver(const std::shared_ptr<AssemblerA>& assemblerA,
                               const std::shared_ptr<AssemblerB>& assemblerB, Eigen::Index nev)
      : PartialGeneralSymEigenSolver(assemblerA->matrix(), assemblerB->matrix(), nev) {
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
   * \return Eigen::VectorXd vector of eigenvalues
   */
  Eigen::VectorXd eigenvalues() const {
    assertCompute();
    return solver_.eigenvalues();
  }

  /**
   * \brief Returns the eigenvectors of the gerneral eigenvalue problem
   *
   * \param _nev optionally specify how many eigenvectors are requested
   * \return auto matrix with the eigevectors as columns
   */
  auto eigenvectors(std::optional<Eigen::Index> _nev = std::nullopt) const {
    assertCompute();
    return solver_.eigenvectors(_nev.value_or(nev_));
  }

  Eigen::Index nev() const { return nev_; }

private:
  Eigen::Index nev_;
  ProductType aOP_;
  CholeskyType bOP_;
  SolverType solver_;

  bool computed_{};

  void assertCompute() const {
    if (not computed_)
      DUNE_THROW(Dune::IOError, "Eigenvalues and -vectors not yet computed, please call compute() first");
  }
};

template <Concepts::FlatAssembler AS1, Concepts::FlatAssembler AS2>
auto makePartialGeneralSymEigenSolver(const std::shared_ptr<AS1>& as1, const std::shared_ptr<AS2> as2) {
  using MatrixType = typename AS1::MatrixType;
  using ScalarType = typename AS1::MatrixType::Scalar;
  using SolverType = PartialGeneralSymEigenSolver<MatrixType>;

  return SolverType{as1, as2};
}

template <Concepts::FlatAssembler AS1, Concepts::FlatAssembler AS2>
PartialGeneralSymEigenSolver(std::shared_ptr<AS1> as1, std::shared_ptr<AS2> as2,
                             int nev) -> PartialGeneralSymEigenSolver<typename AS1::MatrixType>;
} // namespace Ikarus

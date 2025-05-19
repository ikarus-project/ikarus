// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file generalizedeigensolverfactory.hh
 * \brief Factory for classes solving a generalized eigenvalue problem
 */

#pragma once

#include <ikarus/solver/eigenvaluesolver/generalizedeigensolver.hh>

namespace Ikarus {

/**
 * \brief  Factory function to create a GeneralizedSymEigenSolver for a specific backend (Eigen or Spectra) with
 * provided matrices for both quantities.
 *
 * \tparam tag EigenValueSolverType indicating the solver backend.
 * \tparam MATA the type of matrix A
 * \tparam MATB the type of matrix B
 * \return GeneralizedSymEigenSolver The created solver
 */
template <EigenValueSolverType tag, typename MATA, typename MATB>
requires(Concepts::DenseOrSparseEigenMatrix<std::remove_cvref_t<MATA>> &&
         std::same_as<std::remove_cvref_t<MATA>, std::remove_cvref_t<MATB>> &&
         not(tag == EigenValueSolverType::Eigen && Concepts::SparseEigenMatrix<std::remove_cvref_t<MATA>>))
auto makeGeneralizedSymEigenSolver(MATA&& A, MATB&& B) {
  using MatrixType = std::remove_cvref_t<MATA>;
  return GeneralizedSymEigenSolver<tag, MatrixType>{std::forward<MATA>(A), std::forward<MATB>(B)};
}

/**
 * \brief Factory function to create a GeneralizedSymEigenSolver for a specific backend (Eigen or Spectra) with provided
 * assemblers for both quantities.
 *
 * \tparam tag EigenValueSolverType indicating the solver backend.
 * \tparam AssemblerA the type of the assembler for matrix A.
 * \tparam AssemblerB the type of the assembler for matrix B.
 * \param assemblerA assembler for matrix A.
 * \param assemblerB assembler for matrix B.
 * \return GeneralizedSymEigenSolver The created solver
 */
template <EigenValueSolverType tag, Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
requires(std::same_as<typename AssemblerA::MatrixType, typename AssemblerB::MatrixType> &&
         not(tag == EigenValueSolverType::Eigen && Concepts::SparseEigenMatrix<typename AssemblerA::MatrixType>))
auto makeGeneralizedSymEigenSolver(const std::shared_ptr<AssemblerA> assemblerA,
                                   const std::shared_ptr<AssemblerB> assemblerB) {
  makeGeneralizedSymEigenSolver(assemblerA->matrix(), assemblerB->matrix());
}

/**
 * \brief Factory function to create a GeneralizedSymEigenSolver with a provided matrix and an identity matrix.
 * \remark This has some overhead compared to just using the non-generalized versions of the solvers directly, due to
 * the decomposition of the second matrix.
 *
 * \tparam MATA the deduced type of the provided matrix A
 * \param A the provided matrix A
 */
template <EigenValueSolverType tag, typename MATA>
requires(Concepts::DenseOrSparseEigenMatrix<std::remove_cvref_t<MATA>> &&
         not(tag == EigenValueSolverType::Eigen && Concepts::SparseEigenMatrix<std::remove_cvref_t<MATA>>))
auto makeIdentitySymEigenSolver(MATA&& A) {
  using MatrixType = std::remove_cvref_t<MATA>;
  MatrixType I     = MatrixType(A.rows(), A.cols());
  I.setIdentity();

  return GeneralizedSymEigenSolver<tag, MatrixType>{std::forward<MATA>(A), I};
}

/**
 * \brief Factory function to create a GeneralizedSymEigenSolver with a provided matrix from an assembler and an
 * identity matrix.
 *
 * \tparam AssemblerA the type of the assembler for matrix A.
 * \param assemblerA assembler for matrix A.
 */
template <EigenValueSolverType tag, Concepts::FlatAssembler AssemblerA>
requires(Concepts::DenseOrSparseEigenMatrix<typename AssemblerA::MatrixType> &&
         not(tag == EigenValueSolverType::Eigen && Concepts::SparseEigenMatrix<typename AssemblerA::MatrixType>))
auto makeIdentitySymEigenSolver(const std::shared_ptr<AssemblerA>& assemblerA) {
  return makeIdentitySymEigenSolver<tag>(assemblerA->matrix());
}

/**
 * \brief Factory function to create a PartialGeneralizedSymEigenSolver with provided assemblers for both quantities
 *
 * \tparam AssemblerA the type of the assembler for matrix A.
 * \tparam AssemblerB the type of the assembler for matrix B.
 * \param assemblerA assembler for matrix A.
 * \param assemblerB assembler for matrix B.
 * \return PartialGeneralizedSymEigenSolver The created solver
 */
template <Concepts::FlatAssembler AssemblerA, Concepts::FlatAssembler AssemblerB>
requires(std::same_as<typename AssemblerA::MatrixType, typename AssemblerB::MatrixType>)
auto makePartialGeneralizedSymEigenSolver(std::shared_ptr<AssemblerA> assemblerA,
                                          std::shared_ptr<AssemblerB> assemblerB, Eigen::Index nev) {
  return makePartialGeneralizedSymEigenSolver(assemblerA->matrix(), assemblerB->matrix());
}

/**
 * \brief Factory function to create a PartialGeneralizedSymEigenSolver with provided matrices for both quantities
 *
 * \tparam MATA the type of matrix A
 * \tparam MATB the type of matrix B
 * \return PartialGeneralizedSymEigenSolver The created solver
 */
template <typename MATA, typename MATB>
requires(Concepts::DenseOrSparseEigenMatrix<std::remove_cvref_t<MATA>> &&
         std::same_as<std::remove_cvref_t<MATA>, std::remove_cvref_t<MATB>>)
auto makePartialGeneralizedSymEigenSolver(MATA&& A, MATB&& B, Eigen::Index nev) {
  using MatrixType = std::remove_cvref_t<MATA>;
  return PartialGeneralizedSymEigenSolver<MatrixType>{std::forward<MATA>(A), std::forward<MATB>(B), nev};
}

/**
 * \brief Factory function to create a PartialGeneralizedSymEigenSolver with a provided matrix and an identity matrix.
 *
 * \tparam MATA the deduced type of the provided matrix A
 * \param A the provided matrix A
 * \return PartialGeneralizedSymEigenSolver The created solver
 */
template <typename MATA>
requires(Concepts::DenseOrSparseEigenMatrix<std::remove_cvref_t<MATA>>)
auto makePartialIdentitySymEigenSolver(MATA&& A, Eigen::Index nev) {
  using MatrixType = std::remove_cvref_t<MATA>;
  MatrixType I     = MatrixType(A.rows(), A.cols());
  I.setIdentity();

  return PartialGeneralizedSymEigenSolver<MatrixType>{std::forward<MATA>(A), I, nev};
}

/**
 * \brief Factory function to create a PartialGeneralizedSymEigenSolver with a provided matrix from an assembler and an
 * identity matrix.
 *
 * \tparam AssemblerA the type of the assembler for matrix A.
 * \param assemblerA assembler for matrix A.
 */
template <Concepts::FlatAssembler AssemblerA>
auto makePartialIdentitySymEigenSolver(const std::shared_ptr<AssemblerA>& assemblerA, Eigen::Index nev) {
  return makePartialIdentitySymEigenSolver(assemblerA->matrix(), nev);
}

} // namespace Ikarus

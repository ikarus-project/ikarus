//
// Created by alex on 12/22/21.
//

#pragma once
#include <memory>
#include <type_traits>
#include <variant>

#include <Eigen/CholmodSupport>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/UmfPackSupport>
//#include <Eigen/SuperLUSupport>

namespace Ikarus {

  enum class SolverTypeTag {
    si_ConjugateGradient,
    si_LeastSquaresConjugateGradient,
    si_BiCGSTAB,
    sd_SimplicialLLT,
    sd_SimplicialLDLT,
    sd_SparseLU,
    sd_SparseQR,
    sd_CholmodSupernodalLLT,
    sd_UmfPackLU,
    sd_SuperLU,
    d_PartialPivLU,
    d_FullPivLU,
    d_HouseholderQR,
    d_ColPivHouseholderQR,
    d_FullPivHouseholderQR,
    d_CompleteOrthogonalDecomposition,
    d_LLT,
    d_LDLT,
  };

  enum class MatrixTypeTag { Dense, Sparse };

  /** \brief A type-erased solver templated with the scalar type of the linear system */
  template <typename ScalarType = double >
  class ILinearSolver {
  public:
    using SparseMatrixType = Eigen::SparseMatrix<ScalarType>;
    using DenseMatrixType  = Eigen::MatrixX<ScalarType>;
    explicit ILinearSolver(const SolverTypeTag& solverTypeTag) {
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
          throw std::logic_error("Not implemented yet.");
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
          //        case SolverTypeTag::BDCSVD:
          //          solverimpl = std::make_unique<SolverImpl<BDCSVD<DenseMatrixType>>>();
          //          break;
          //        case SolverTypeTag::JacobiSVD:
          //          solverimpl = std::make_unique<SolverImpl<JacobiSVD<DenseMatrixType>>>();
          break;
        case SolverTypeTag::d_LDLT:
          solverimpl = std::make_unique<SolverImpl<LDLT<DenseMatrixType>>>();
          break;
      }
    }

    ~ILinearSolver()       = default;
    ILinearSolver& operator=(const ILinearSolver& other) {
      ILinearSolver tmp(other);
      std::swap(solverimpl, tmp.solverimpl);
      return *this;
    }

    ILinearSolver(ILinearSolver&&) noexcept = default;
    ILinearSolver& operator=(ILinearSolver&&) noexcept = default;

  private:
    struct SolverBase {
      virtual ~SolverBase() = default;
      //    [[nodiscard]] virtual std::unique_ptr<SolverBase> clone() const                                  = 0;
      virtual void analyzePattern(const DenseMatrixType&) const {};
      virtual void analyzePattern(const SparseMatrixType&)                              = 0;
      virtual void factorize(const DenseMatrixType&)                                    = 0;
      virtual void factorize(const SparseMatrixType&)                                   = 0;
      virtual void compute(const SparseMatrixType&)                                     = 0;
      virtual void compute(const DenseMatrixType&)                                      = 0;
      virtual void solve(Eigen::VectorX<ScalarType>&x, const Eigen::VectorX<ScalarType>&) const = 0;
    };

    template <typename Solver>
    struct SolverImpl : public SolverBase {

      void analyzePattern(const SparseMatrixType& A) override {
        if constexpr (requires(Solver sol) { sol.analyzePattern(A); }) solver.analyzePattern(A);
      }

      void factorize(const SparseMatrixType& A) override {
        if constexpr (requires(Solver sol) { sol.factorize(A); })
          solver.factorize(A);
      }

      // Dense Solvers do not have a factorize method therefore
      // our interface we just call compute for dense matrices
      void factorize(const DenseMatrixType& A) override {
        if constexpr (requires(Solver sol) { sol.compute(A); } && std::is_base_of_v<Eigen::SolverBase<Solver>, Solver>)
          solver.compute(A);
      }
      void compute(const SparseMatrixType& A) {
        //        std::cout << "solver.computeSparse(A)" << std::endl;
        if constexpr (std::is_base_of_v<Eigen::SparseSolverBase<Solver>, Solver>)
          solver.compute(A);
        else
          throw std::logic_error("This solver does not support solving with sparse matrices.");
      }

      void compute(const DenseMatrixType& A) {
        if constexpr (std::is_base_of_v<Eigen::SolverBase<Solver>, Solver>)
          solver.compute(A);
        else
          throw std::logic_error("This solver does not support solving with dense matrices.");
      }

      void solve(Eigen::VectorX<ScalarType>&x, const Eigen::VectorX<ScalarType>& b) const override {
        x = solver.solve(b);
      }

      Solver solver;
    };

    std::unique_ptr<SolverBase> solverimpl;

  public:
    template <typename MatrixType>
    requires std::is_same_v<MatrixType, DenseMatrixType> || std::is_same_v<MatrixType, SparseMatrixType>
    inline ILinearSolver& compute(const MatrixType& A) {
      //            std::cout << "compute(A)" << std::endl;
      //            std::cout <<"r,c: "<< A.rows()<<" "<<A.cols() << std::endl;
      //            std::cout << A << std::endl;
      solverimpl->compute(A);
      return *this;
    }
    template <typename MatrixType>
    requires std::is_same_v<MatrixType, DenseMatrixType> || std::is_same_v<MatrixType, SparseMatrixType>
    inline void analyzePattern(const MatrixType& A) {
      //      std::cout << "analyzePattern(A)" << std::endl;
      solverimpl->analyzePattern(A);
    }

    template <typename MatrixType>
    requires std::is_same_v<MatrixType, DenseMatrixType> || std::is_same_v<MatrixType, SparseMatrixType>
    inline void factorize(const MatrixType& A) {
      //      std::cout << "factorize(A)" << std::endl;
      solverimpl->factorize(A);
    }

    void solve(Eigen::VectorX<ScalarType>&x, const Eigen::VectorX<ScalarType>& b) {
      solverimpl->solve(x,b);
    }
  };

}  // namespace Ikarus
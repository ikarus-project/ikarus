//
// Created by alex on 12/22/21.
//

#pragma once
#include <type_traits>
#include <memory>
#include <variant>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
#include<Eigen/SparseQR>
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
//#include <Eigen/SuperLUSupport>


namespace Ikarus
{




enum class SolverTypeTag { CG, LeastSquareCG,BiCGSTAB,SimplicialLLT,SimplicialLDLT,SparseLU,SparseQR,CholmodSupernodalLLT,UmfPackLU,SuperLU,PartialPivLU,FullPivLU,HouseholderQR,ColPivHouseholderQR,FullPivHouseholderQR,CompleteOrthogonalDecomposition,LLT,LDLT,BDCSVD,JacobiSVD};
  template<SolverTypeTag solverTypeTag, typename ScalarType=double, typename... SolverOptions>
  class LinearSolveFactory{

   static auto createSolver()
    {
     if constexpr (solverTypeTag==SolverTypeTag::CG) {
       if constexpr (sizeof...(SolverOptions) == 0)
         return Eigen::ConjugateGradient<Eigen::SparseMatrix<ScalarType>, Eigen::Lower | Eigen::Upper>();
     } else if constexpr (solverTypeTag==SolverTypeTag::SimplicialLDLT)
        return  Eigen::SimplicialLDLT<Eigen::SparseMatrix<ScalarType>>();
   }
  };





template<typename SolverType, typename ScalarType=double>
class DenseLinearSolver {
public:

  using MatrixType =Eigen::MatrixX<ScalarType>;

  void  analyzePattern(const MatrixType& mat){  };



};

enum class MatrixTypeTag {Dense,Sparse};

/** \brief A type-erased solver templated with the matix type */
template<typename ScalarType>
class ILinearSolver {
public:
  template <typename SolverType>
  explicit ILinearSolver(const SolverType &fe) : feimpl{std::make_unique<SolverImpl<SolverType>>(fe)} {
  }

  ~ILinearSolver() = default;
  ILinearSolver(const ILinearSolver &other) : feimpl{other.feimpl->clone()} {}
  ILinearSolver &operator=(const ILinearSolver &other) {
    ILinearSolver tmp(other);
    std::swap(feimpl, tmp.feimpl);
    return *this;
  }

  ILinearSolver(ILinearSolver &&) noexcept = default;
  ILinearSolver &operator=(ILinearSolver &&) noexcept = default;

  using SparseMatrixType = Eigen::SparseMatrix<ScalarType>;
  using DenseMatrixType = Eigen::MatrixX<ScalarType>;

private:
  struct SolverBase {
    virtual ~SolverBase()                            = default;
    [[nodiscard]] virtual std::unique_ptr<SolverBase> clone() const                                  = 0;
    virtual void analyzePattern(const DenseMatrixType&) const { };
    virtual void analyzePattern(const SparseMatrixType&)    = 0;
    virtual void factorize(const DenseMatrixType&) const { };
    virtual void factorize(const SparseMatrixType&) const    = 0;
    virtual Eigen::VectorX<ScalarType> solve(const Eigen::VectorX<ScalarType>&) const    = 0;
  };

  template <typename Solver>
  struct SolverImpl : public SolverBase {
    explicit SolverImpl(Solver solverarg) : solver{solverarg} {};
    void analyzePattern(const SparseMatrixType& A)   override  {solver.analyzePattern(A);}
    void factorize(const SparseMatrixType& A)   override  {solver.factorize(A);}
    Eigen::VectorX<ScalarType> solve(const Eigen::VectorX<ScalarType>& b) const   override  {solver.solve(b);}

    [[nodiscard]] std::unique_ptr<SolverBase> clone() const final { return std::make_unique<SolverImpl>(*this); }
    Solver solver;
  };

  std::unique_ptr<SolverBase> feimpl;

public:
  template <typename MatrixType> requires std::is_same_v<MatrixType,DenseMatrixType> ||  std::is_same_v<MatrixType,SparseMatrixType>
  auto inline compute(const MatrixType& A){
    if constexpr (std::is_same_v<MatrixType,SparseMatrixType> ) {
      feimpl->analyzePattern(A);
      feimpl->factorize(A);
    }

    return *this;
  }

  [[nodiscard]] Eigen::VectorX<ScalarType> solve(const Eigen::VectorX<ScalarType>& b)
  {
    return feimpl->solve(b);
  }

};





template<typename SolverType, typename ScalarType=double>
class SparseLinearSolver {
public:

  using MatrixType =Eigen::MatrixX<ScalarType>;

  void  analyzePattern(const MatrixType& mat){
    analyzePatternImpl(mat);
  };
};
template<typename SolverType, typename ScalarType=double>
using LinearSolver = std::variant<SparseLinearSolver<SolverType,ScalarType>,DenseLinearSolver<SolverType,ScalarType>>;


}
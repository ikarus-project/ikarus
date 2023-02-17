// Original File: https://gitlab.com/libeigen/eigen/-/blob/master/Eigen/src/IterativeLinearSolvers/ConjugateGradient.h
// SPDX-FileCopyrightText: 2011-2014 Copyright (C) Gael Guennebaud <gael.guennebaud@inria.fr>
// SPDX-License-Identifier: MPL-2.0
// Modifications:
// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

// this file is modified from the default conjugate gradient method of eigen
#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen {
  enum class TCGStopReason : int {
    negativeCurvature,
    exceededTrustRegion,
    reachedTargetResidualKappaLinear,
    reachedTargetResidualThetaSuperLinear,
    maximumInnerIterations,
    modelIncreased
  };
  template <typename Scalar>
  struct TCGInfo {
    TCGStopReason stop_tCG = TCGStopReason::maximumInnerIterations;
    Scalar Delta           = 100000;
    Scalar kappa           = 0.1;
    Scalar theta           = 1.0;
    int mininner           = 1;
    int maxinner           = 1000;
    int numInnerIter       = 0;

    void initRuntimeOptions(int _num_dof_solve) {
      maxinner = _num_dof_solve * 2;
      Delta    = 100000;  // typical distance of manifold
      //            Delta0 = Delta_bar/8.0;
    }
  };
  namespace internal {

    /** \internal Low-level conjugate gradient algorithm
     * \param mat The matrix A
     * \param rhs The right hand side vector b
     * \param x On input and initial solution, on output the computed solution.
     * \param precond A preconditioner being able to efficiently solve for an
     *                approximation of Ax=b (regardless of b)
     * \param iters On input the max number of iteration, on output the number of performed iterations.
     * \param tol_error On input the tolerance error, on output an estimation of the relative error.
     */
    template <typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
    void truncated_conjugate_gradient(const MatrixType& mat, const Rhs& rhs, Dest& x, const Preconditioner& precond,
                                      Eigen::Index& iters, typename Dest::RealScalar& tol_error,
                                      TCGInfo<typename Dest::RealScalar>& _info) {
      using std::abs;
      using std::sqrt;
      typedef typename Dest::RealScalar RealScalar;
      typedef typename Dest::Scalar Scalar;
      typedef Matrix<Scalar, Dynamic, 1> VectorType;

      RealScalar tol = tol_error;
      Index maxIters = iters;

      Index n = mat.cols();

      VectorType residual = rhs - mat * x;  // initial residual

      RealScalar rhsNorm2 = rhs.norm();
      if (rhsNorm2 == 0) {
        x.setZero();
        iters     = 0;
        tol_error = 0;
        return;
      }
      const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();
      RealScalar threshold            = numext::maxi(tol * tol * rhsNorm2 * rhsNorm2, considerAsZero);
      RealScalar residualNorm2        = residual.norm();
      if (residualNorm2 * residualNorm2 < threshold) {
        iters     = 0;
        tol_error = (residualNorm2 / rhsNorm2);
        return;
      }

      double e_Pd     = 0.0;
      double e_Pe_new = 0.0;
      double e_Pe     = x.squaredNorm();
      double d_Pd;
      double d_Hd;
      VectorType p(n);
      p = precond.solve(residual);  // initial search direction
                                    // bool coutflag=true;
      _info.stop_tCG = TCGStopReason::maximumInnerIterations;
      VectorType z(n), tmp(n);
      RealScalar absNew = numext::real(residual.dot(p));  // the square of the absolute value of r scaled by invM
      d_Pd              = absNew;
      Index i           = 1;
      while (i < maxIters) {
        tmp.noalias() = mat * p;  // the bottleneck of the algorithm
        d_Hd          = p.dot(tmp);
        Scalar alpha  = absNew / d_Hd;  // the amount we travel on dir

        e_Pe_new = e_Pe + 2.0 * alpha * e_Pd + alpha * alpha * d_Pd;

        if (d_Hd <= 0 || e_Pe_new >= _info.Delta * _info.Delta)  // negative curvature or execdet trustregion
        {
          double tau = (-e_Pd + sqrt(e_Pd * e_Pd + d_Pd * (_info.Delta * _info.Delta - e_Pe))) / d_Pd;

          x += tau * p;
          if (d_Hd <= 0)
            _info.stop_tCG = TCGStopReason::negativeCurvature;
          else
            _info.stop_tCG = TCGStopReason::exceededTrustRegion;

          break;
        }
        e_Pe = e_Pe_new;
        x += alpha * p;           // update solution
        residual -= alpha * tmp;  // update residual

        residualNorm2 = residual.norm();

        if (i >= _info.mininner
            && residualNorm2 <= rhsNorm2 * std::min(rhsNorm2, _info.kappa))  // missing pow(rhsNorm2,_info.theta
        {
          // Residual is small enough to quit
          if (_info.kappa < rhsNorm2)
            _info.stop_tCG = TCGStopReason::reachedTargetResidualKappaLinear;
          else
            _info.stop_tCG = TCGStopReason::reachedTargetResidualThetaSuperLinear;
          break;
        }
        if (residualNorm2 < threshold) break;

        z = precond.solve(residual);  // approximately solve for "A z = residual"

        RealScalar absOld = absNew;
        absNew            = numext::real(residual.dot(z));  // update the absolute value of r
        RealScalar beta = absNew / absOld;  // calculate the Gram-Schmidt value used to create the new search direction

        e_Pd = beta * (e_Pd + alpha * d_Pd);
        d_Pd = absNew + beta * beta * d_Pd;

        p = z + beta * p;  // update search direction
        i++;
      }
      tol_error          = (residualNorm2 / rhsNorm2);
      iters              = i;
      _info.numInnerIter = i;
    }

  }  // namespace internal

  template <typename _MatrixType, int _UpLo = Lower,
            typename _Preconditioner = DiagonalPreconditioner<typename _MatrixType::Scalar> >
  class TruncatedConjugateGradient;

  namespace internal {

    template <typename _MatrixType, int _UpLo, typename _Preconditioner>
    struct traits<TruncatedConjugateGradient<_MatrixType, _UpLo, _Preconditioner> > {
      typedef _MatrixType MatrixType;
      typedef _Preconditioner Preconditioner;
    };

  }  // namespace internal

  /** \ingroup IterativeLinearSolvers_Module
      * \brief A conjugate gradient solver for sparse (or dense) self-adjoint problems
      *
      * This class allows to solve for A.x = b linear problems using an iterative conjugate gradient algorithm.
      * The matrix A must be selfadjoint. The matrix A and the vectors x and b can be either dense or sparse.
      *
      * \tparam _MatrixType the type of the matrix A, can be a dense or a sparse matrix.
      * \tparam _UpLo the triangular part that will be used for the computations. It can be Lower,
      *               \c Upper, or \c Lower|Upper in which the full matrix entries will be considered.
      *               Default is \c Lower, best performance is \c Lower|Upper.
      * \tparam _Preconditioner the type of the preconditioner. Default is DiagonalPreconditioner
      *
      * \implsparsesolverconcept
      *
      * The maximal number of iterations and tolerance value can be controlled via the setMaxIterations()
      * and setTolerance() methods. The defaults are the size of the problem for the maximal number of iterations
      * and NumTraits<Scalar>::epsilon() for the tolerance.
      *
      * The tolerance corresponds to the relative residual error: |Ax-b|/|b|
      *
      * \b Performance: Even though the default value of \c _UpLo is \c Lower, significantly higher performance is
      * achieved when using a complete matrix and \b Lower|Upper as the \a _UpLo template parameter. Moreover, in this
      * case multi-threading can be exploited if the user code is compiled with OpenMP enabled.
      * See \ref TopicMultiThreading for details.
      *
      * This class can be used as the direct solver classes. Here is a typical usage example:
        \code
        int n = 10000;
        VectorXd x(n), b(n);
        SparseMatrix<double> A(n,n);
        // fill A and b
        ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
        cg.compute(A);
        x = cg.solve(b);
        std::cout << "#iterations:     " << cg.iterations() << std::endl;
        std::cout << "estimated error: " << cg.error()      << std::endl;
        // update b, and solve again
        x = cg.solve(b);
        \endcode
      *
      * By default the iterations start with x=0 as an initial guess of the solution.
      * One can control the start using the solveWithGuess() method.
      *
      * ConjugateGradient can also be used in a matrix-free context, see the following \link MatrixfreeSolverExample
     example \endlink.
      *
      * \sa class LeastSquaresConjugateGradient, class SimplicialCholesky, DiagonalPreconditioner,
     IdentityPreconditioner
      */
  template <typename _MatrixType, int _UpLo, typename _Preconditioner>
  class TruncatedConjugateGradient
      : public IterativeSolverBase<TruncatedConjugateGradient<_MatrixType, _UpLo, _Preconditioner> > {
  public:
    typedef IterativeSolverBase<TruncatedConjugateGradient> Base;
    TruncatedConjugateGradient(TruncatedConjugateGradient&& other) noexcept
        : Base(std::move(other)), algInfo{other.algInfo} {};

  private:
    using Base::m_error;
    using Base::m_info;
    using Base::m_isInitialized;
    using Base::m_iterations;
    using Base::matrix;
    mutable TCGInfo<typename _MatrixType::RealScalar> algInfo;

  public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef _Preconditioner Preconditioner;

    enum { UpLo = _UpLo };

  public:
    TCGInfo<typename _MatrixType::RealScalar> getInfo() { return algInfo; }

    void setInfo(TCGInfo<typename _MatrixType::RealScalar> _alginfo) { this->algInfo = _alginfo; }
    /** Default constructor. */
    TruncatedConjugateGradient() : Base() {}

    /** Initialize the solver with matrix \a A for further \c Ax=b solving.
     *
     * This constructor is a shortcut for the default constructor followed
     * by a call to compute().
     *
     * \warning this class stores a reference to the matrix A as well as some
     * precomputed values that depend on it. Therefore, if \a A is changed
     * this class becomes invalid. Call compute() to update it with the new
     * matrix A, or modify a copy of A.
     */
    template <typename MatrixDerived>
    explicit TruncatedConjugateGradient(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}

    ~TruncatedConjugateGradient() {}

    /** \internal */
    template <typename Rhs, typename Dest>
    void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const {
      typedef typename Base::MatrixWrapper MatrixWrapperBase;
      typedef typename Base::ActualMatrixType ActualMatrixTypeBase;
      enum {
        TransposeInput = (!MatrixWrapperBase::MatrixFree) && (UpLo == (Lower | Upper)) && (!MatrixType::IsRowMajor)
                         && (!NumTraits<Scalar>::IsComplex)
      };
      typedef typename internal::conditional<TransposeInput, Transpose<const ActualMatrixTypeBase>,
                                             ActualMatrixTypeBase const&>::type RowMajorWrapper;
      EIGEN_STATIC_ASSERT(EIGEN_IMPLIES(MatrixWrapperBase::MatrixFree, UpLo == (Lower | Upper)),
                          MATRIX_FREE_CONJUGATE_GRADIENT_IS_COMPATIBLE_WITH_UPPER_UNION_LOWER_MODE_ONLY);
      typedef typename internal::conditional<
          UpLo == (Lower | Upper), RowMajorWrapper,
          typename MatrixWrapperBase::template ConstSelfAdjointViewReturnType<UpLo>::Type>::type SelfAdjointWrapper;
      //      TCGInfo<typename Dest::RealScalar> dummyInfo;
      m_iterations = Base::maxIterations();
      m_error      = Base::m_tolerance;

      RowMajorWrapper row_mat(matrix());
      internal::truncated_conjugate_gradient(SelfAdjointWrapper(row_mat), b, x, Base::m_preconditioner, m_iterations,
                                             m_error, algInfo);
      m_info = m_error <= Base::m_tolerance ? Success : NoConvergence;

      //      setInfo =dummyInfo;
    }

  protected:
  };

}  // end namespace Eigen

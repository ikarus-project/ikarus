// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <Spectra/SymEigsSolver.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ikarus/linearalgebra/nonLinearoperator.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/fancyactivationfunctions.hh>

namespace Ikarus {

  template <typename ActivationFunction = Ikarus::AdaptiveStepSizing::LinearPiecewiseFunction,
            typename ScalarType         = double>
  struct DefaultAdaptiveStepSizing {
    template <typename NonLinearOperator>
    void operator()(const NonLinearSolverInformation& solverInfo, const NonLinearOperator&,
                    SubsidiaryArgs& subsidiaryArgs, const ActivationFunction&) {
      if (targetIterations != 0 && subsidiaryArgs.actStep > 0) {
        ScalarType actIte       = solverInfo.iterations;
        subsidiaryArgs.stepSize = sqrt(static_cast<ScalarType>(targetIterations) / actIte) * subsidiaryArgs.stepSize;
      }
    }
    int targetIterations{0};
  };

  template <typename ActivationFunction = Ikarus::AdaptiveStepSizing::LinearPiecewiseFunction,
            typename ScalarType         = double>
  struct FancyAdaptiveStepSizing {
    template <typename NonLinearOperator>
    void operator()(const NonLinearSolverInformation& solverInfo, NonLinearOperator& nonLinearOperator,
                    SubsidiaryArgs& subsidiaryArgs, const ActivationFunction& activationFunction) const {
      const auto& K = nonLinearOperator.derivative();
      // const  auto& Kred = K.getReducedMatrixImpl();
      auto act_step = subsidiaryArgs.actStep;

      /// calculate dimension of K (stiffness matrix)
      auto K_size = K.rows();

      /// Construct matrix operation object using the wrapper class DenseSymMatProd
      Spectra::DenseSymMatProd<double> op(K);

      /// Construct eigen solver object, requesting the largest three eigenvalues
      /// Refer: https://spectralib.org/doc/classspectra_1_1symeigssolver#a469d37e08be076e198c9e383236c5056
      auto num_eigs  = 3;
      auto eigs_conv = 2 * num_eigs;
      assert((eigs_conv > num_eigs) and (eigs_conv <= K_size)
             && "eigs_conv should be chosen such that (num_eigs < eigs_conv <= K_size)");

      Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, num_eigs, eigs_conv);

      /// Initialize and compute
      eigs.init();
      // eigs.compute(Spectra::SortRule::SmallestAlge);
      eigs.compute(Spectra::SortRule::LargestAlge);  // TODO: Change it to SmallestAlge
      Eigen::VectorXd evalues(num_eigs);
      if (eigs.info() == Spectra::CompInfo::Successful) evalues = eigs.eigenvalues();

      std::cout << "Eigenvalues found:\n" << evalues << std::endl;

      std::vector<double> evstd;
      for (auto i = 0U; i < evalues.size(); ++i)
        if (std::abs(evalues[i] - 1.0)
            < 1e-8)  // TODO: Fix it such that it works for both reduced and full matrices correctly
          evstd.push_back(evalues[i] * 1e10);
        else
          evstd.push_back(evalues[i]);

      std::sort(evstd.begin(), evstd.end(), std::less<>());
      for (auto i = 0U; i < evalues.size(); ++i)
        evalues[i] = evstd[i];

      std::cout << "Eigenvalues sorted:\n" << evalues << std::endl;

      /// store eigenvalues of each increment in a matrix
      subsidiaryArgs.eigenValues.col(act_step) = evalues;
      auto& ev                                 = subsidiaryArgs.eigenValues;

      if (act_step > 0) {
        /// calculation r_soll with piecewise linear function
        double lambdaHat = ev(0, act_step) / ev(0, 0);
        double r_soll    = activationFunction(lambdaHat);
        /// calculation scale factor
        double factor_x{0.0};
        factor_x = r_soll * ev(0, 0) / (ev(0, act_step - 1) - ev(0, act_step));
        /// scaled step size
        subsidiaryArgs.stepSize = factor_x * subsidiaryArgs.stepSize;
      }
    }
  };
}  // namespace Ikarus

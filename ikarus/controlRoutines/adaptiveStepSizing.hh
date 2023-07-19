// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <Spectra/SymEigsSolver.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "ikarus/solver/nonLinearSolver/newtonRaphson.hh"
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphsonWithScalarSubsidiaryFunction.hh>

namespace Ikarus {

  // TODO move these structs to seperate header
  //  sync solverinformations ,NR,NRwithSubArgs,TrustRegion
  //  improve output
  template <typename ScalarType = double>
  struct DefaultAdaptiveStepSizing {
    template <typename NonLinearOperator>
    void operator()(const SolverInfos& solverInfo, const NonLinearOperator&, SubsidiaryArgs& subsidiaryArgs) {
      if (targetIterations != 0 && subsidiaryArgs.actStep > 0) {
        ScalarType actIte       = solverInfo.iterations;
        subsidiaryArgs.stepSize = sqrt(static_cast<ScalarType>(targetIterations) / actIte) * subsidiaryArgs.stepSize;
      }
    }
    int targetIterations{0};
  };

  template <typename ScalarType = double>
  struct FancyAdaptiveStepSizing {
    template <typename NonLinearOperator>
    void operator()(const SolverInfos& solverInfo, NonLinearOperator& nonLinearOperator,
                    SubsidiaryArgs& subsidiaryArgs) const {
      const auto& K = nonLinearOperator.derivative();
      auto act_step = subsidiaryArgs.actStep;

      /// Construct matrix operation object using the wrapper class DenseSymMatProd
      Spectra::DenseSymMatProd<double> op(K);

      /// Construct eigen solver object, requesting the largest three eigenvalues
      auto num_eigs  = 3;  // TODO: Generalize it
      auto eigs_conv = 2 * num_eigs;
      Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, num_eigs, eigs_conv);

      /// Initialize and compute
      eigs.init();
      eigs.compute(Spectra::SortRule::LargestAlge);  // TODO: Change it to SmallestAlge
      Eigen::VectorXd evalues(num_eigs);
      if (eigs.info() == Spectra::CompInfo::Successful) evalues = eigs.eigenvalues();

      std::vector<double> evstd;
      for (auto i = 0U; i < evalues.size(); ++i)
        if (std::abs(evalues[i] - 1.0)
            < 1e-8)  // TODO: Fix it such that it works for both reduced and full matrices correctly
          evstd.push_back(evalues[i] * 1e10);
        else
          evstd.push_back(evalues[i]);

      std::sort(evstd.begin(), evstd.end(), std::less<double>());
      for (auto i = 0U; i < evalues.size(); ++i)
        evalues[i] = evstd[i];

      /// store eigenvalues of each increment in a matrix
      subsidiaryArgs.eigenValues.col(act_step) = evalues;
      auto& ev                                 = subsidiaryArgs.eigenValues;

      if (act_step > 0) {
        /// calculation r_soll with piecewise linear function
        double r_soll     = 0.0;
        double r_0        = 0.01;
        double r_i        = 0.05;
        double r_1        = 0.2;
        double lambda_i   = 0.1;
        double lambda_hat = ev(0, act_step) / ev(0, 0);
        if (lambda_hat <= lambda_i) {
          r_soll = r_0 + (r_i - r_0) / lambda_i * lambda_hat;
        } else {
          r_soll = (lambda_i * r_1 - r_i) / (lambda_i - 1) + (r_1 - r_i) / (1 - lambda_i) * lambda_hat;
        }

        /// calculation scale factor
        double factor_x{0.0};
        factor_x = r_soll * ev(0, 0) / (ev(0, act_step - 1) - ev(0, act_step));
        /// scaled step size
        subsidiaryArgs.stepSize = factor_x * subsidiaryArgs.stepSize;
      }
    }
    int targetIterations{0};
  };
}  // namespace Ikarus

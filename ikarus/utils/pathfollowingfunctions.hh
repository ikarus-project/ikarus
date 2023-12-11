// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <math.h>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <Eigen/Core>

#include <ikarus/solver/linearsolver/linearsolver.hh>

namespace Ikarus {

  struct SubsidiaryArgs {
    double stepSize;
    Eigen::VectorX<double> DD;
    double Dlambda{};
    double f{};
    Eigen::VectorX<double> dfdDD;
    double dfdDlambda{};
  };

  /// Arc Length Control Method
  struct StandardArcLength {
    void evaluateSubsidiaryFunction(SubsidiaryArgs& args) const {
      if (psi) {
        const auto root = sqrt(args.DD.squaredNorm() + psi.value() * psi.value() * args.Dlambda * args.Dlambda);
        args.f          = root - args.stepSize;
        args.dfdDD      = args.DD / root;
        args.dfdDlambda = (psi.value() * psi.value() * args.Dlambda) / root;
      } else
        DUNE_THROW(Dune::InvalidStateException,
                   "You have to call initialPrediction first. Otherwise psi is not defined");
    }

    template <typename NonLinearOperator>
    void initialPrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args) {
      auto linearSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);  // for the linear predictor step

      nonLinearOperator.lastParameter() = 1.0;  // lambda =1.0

      nonLinearOperator.template update<0>();
      const auto& R = nonLinearOperator.value();
      const auto& K = nonLinearOperator.derivative();

      linearSolver.factorize(K);
      linearSolver.solve(args.DD, -R);

      const auto DD2 = args.DD.squaredNorm();

      psi    = sqrt(DD2);
      auto s = sqrt(psi.value() * psi.value() + DD2);

      args.DD      = args.DD * args.stepSize / s;
      args.Dlambda = args.stepSize / s;

      nonLinearOperator.firstParameter() = args.DD;
      nonLinearOperator.lastParameter()  = args.Dlambda;
    }

    template <typename NonLinearOperator>
    void intermediatePrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args) {
      nonLinearOperator.firstParameter() += args.DD;
      nonLinearOperator.lastParameter() += args.Dlambda;
    }

    std::string name = "Arc length";

  private:
    std::optional<double> psi;
  };

  /// Load Control Method
  struct LoadControlWithSubsidiaryFunction {
    void evaluateSubsidiaryFunction(SubsidiaryArgs& args) const {
      args.f = args.Dlambda - args.stepSize;
      args.dfdDD.setZero();
      args.dfdDlambda = 1.0;
    }

    template <typename NonLinearOperator>
    void initialPrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args) {
      args.Dlambda                      = args.stepSize;
      nonLinearOperator.lastParameter() = args.Dlambda;
    }

    template <typename NonLinearOperator>
    void intermediatePrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args) {
      nonLinearOperator.lastParameter() += args.Dlambda;
    }
    std::string name = "Load control";
  };

  /// Displacement Control Method
  struct DisplacementControl {
    explicit DisplacementControl(std::vector<int> p_controlledIndices)
        : controlledIndices{std::move(p_controlledIndices)} {}

    void evaluateSubsidiaryFunction(SubsidiaryArgs& args) const {
      auto controlledDOFsNorm = args.DD(controlledIndices).norm();
      args.f                  = controlledDOFsNorm - args.stepSize;
      args.dfdDlambda         = 0.0;
      args.dfdDD.setZero();
      args.dfdDD(controlledIndices) = args.DD(controlledIndices) / controlledDOFsNorm;
    }

    template <typename NonLinearOperator>
    void initialPrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args) {
      args.DD(controlledIndices).array() = args.stepSize;
      nonLinearOperator.firstParameter() = args.DD;
    }

    template <typename NonLinearOperator>
    void intermediatePrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args) {
      nonLinearOperator.firstParameter() += args.DD;
    }

    std::string name = "Displacement control";

  private:
    std::vector<int> controlledIndices;
  };

}  // namespace Ikarus

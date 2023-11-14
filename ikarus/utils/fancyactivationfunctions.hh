// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <memory>
#include <optional>

#include <dune/common/float_cmp.hh>

#include <ikarus/solver/linearsolver/linearsolver.hh>

namespace Ikarus::AdaptiveStepSizing {

  struct ActivationFunctionArgs {
    double r_0{0.01};
    double r_i{0.05};
    double r_1{0.2};
    double lambda_i{0.1};
  };

  /// Constant Function
  struct ConstantFunction {
    ConstantFunction() = default;
    explicit ConstantFunction(ActivationFunctionArgs& p_args) { args = p_args; }
    double operator()(double) const { return args.r_0; }
    std::string name = "Constant Function";
    ActivationFunctionArgs args;
  };

  /// Linear Piecewise Function
  struct LinearPiecewiseFunction {
    LinearPiecewiseFunction() = default;
    explicit LinearPiecewiseFunction(ActivationFunctionArgs& p_args) { args = p_args; }
    double operator()(double lambda_hat) const {
      double r_soll = 0.0;
      if (lambda_hat <= args.lambda_i) {
        r_soll = args.r_0 + (args.r_i - args.r_0) / args.lambda_i * lambda_hat;
      } else {
        r_soll = (args.lambda_i * args.r_1 - args.r_i) / (args.lambda_i - 1)
                 + (args.r_1 - args.r_i) / (1 - args.lambda_i) * lambda_hat;
      }
      return r_soll;
    }
    std::string name = "Linear Piecewise Function";
    ActivationFunctionArgs args;
  };

}  // namespace Ikarus::AdaptiveStepSizing

// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testCommon.hh"

#include <variant>

#include <ikarus/finiteElements/mechanics/enhancedAssumedStrains.hh>

template <typename DisplacementBasedElement>
struct ElementTest<Ikarus::EnhancedAssumedStrains<DisplacementBasedElement>> {
  [[nodiscard]] static auto test() {
    auto easFunctor = [](auto& nonLinOp, auto& fe, auto& req) {
      const auto& localView = fe.localView();
      const auto& element   = localView.element();
      constexpr int gridDim = std::remove_cvref_t<decltype(element)>::dimension;

      Dune::TestSuite t("EAS specific test");

      auto subOp = nonLinOp.template subOperator<1, 2>();
      std::array<int, (gridDim == 2) ? 4 : 3> easParameters;
      if constexpr (gridDim == 2)
        easParameters = {0, 4, 5, 7};
      else
        easParameters = {0, 9, 21};

      Eigen::VectorXd oldEigenValues, newEigenValues;

      for (auto& numberOfEASParameter : easParameters) {
        fe.setEASType(numberOfEASParameter);
        subOp.updateAll();
        auto messageIfFailed = "The numbers of EAS parameters are " + std::to_string(numberOfEASParameter) + ".";
        if (numberOfEASParameter == 0) {
          nonLinOp.updateAll();
          t.subTest(checkGradientOfElement(nonLinOp, messageIfFailed));
          t.subTest(checkHessianOfElement(nonLinOp, messageIfFailed));
        }
        t.subTest(checkJacobianOfElement(subOp, messageIfFailed));
        t.subTest(checkFEByAutoDiff(nonLinOp, fe, req, messageIfFailed));

        auto stiffnessMatrix = subOp.derivative();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(stiffnessMatrix);
        newEigenValues = es.eigenvalues();

        t.check((newEigenValues.array() < 1e-6 * newEigenValues.norm()).count() == 3 * gridDim - 3,
                "Number of eigenvalues")
            << "We always should have " << 3 * gridDim - 3 << " zero eigenvalues, for " << 3 * gridDim - 3
            << " rigid body motions"
            << " but we counted " << (newEigenValues.array() < 1e-6 * newEigenValues.norm()).count()
            << "\nEigenValues: \n"
            << newEigenValues.transpose() << "The tolerance is " << 1e-6 * newEigenValues.norm()
            << "The number of EAS parameters is" << numberOfEASParameter << std::endl;
        if (numberOfEASParameter > 0 and numberOfEASParameter != 5)          // Q1E4 and Q1E5 are the same
          t.check((newEigenValues.array() <= oldEigenValues.array()).sum())  // Q1E4 and Q1E7 are the same in the case
                                                                             // of undistorted element
              << "More EAS parameter mean that the stiffness gets reduced. EAS parameter: " << numberOfEASParameter
              << "\noldEigenValues: \n"
              << oldEigenValues.transpose() << "\nnewEigenValues: \n"
              << newEigenValues.transpose() << std::endl;
        oldEigenValues = newEigenValues;

        const int order = 2 * (localView.tree().child(0).finiteElement().localBasis().order());
        auto rule       = Dune::QuadratureRules<double, gridDim>::rule(localView.element().type(), order);

        if (numberOfEASParameter > 0) {
          auto easVariantCopy    = fe.easVariant();  // This only test if the variant has a copy assignment operator
          const auto& easVariant = fe.easVariant();
          std::visit(
              [&]<typename EAS>(const EAS& easFunction) {
                if constexpr (not std::is_same_v<std::monostate, EAS>) {
                  typename EAS::MType MIntegrated;
                  MIntegrated.setZero();
                  for (const auto& gp : rule) {
                    const auto M = easFunction.calcM(gp.position());

                    const double detJ = element.geometry().integrationElement(gp.position());
                    MIntegrated += M * detJ * gp.weight();
                  }
                  t.check(MIntegrated.isZero())
                      << "Orthogonality condition check: The M matrix of the EAS method should be "
                         "zero, integrated over the domain.";
                }
              },
              easVariant);
          try {
            auto requirements = Ikarus::FErequirements();
            fe.calculateScalar(requirements);
            t.check(false) << "fe.calculateScalar should have failed for numberOfEASParameter > 0";
          } catch (const Dune::NotImplemented&) {
          }
        }

        try {
          fe.setEASType(100);
          t.check(false) << "fe.setEASType(100) should have failed";
        } catch (const Dune::NotImplemented&) {
        }
      };
      return t;
    };
    return easFunctor;
  }
};

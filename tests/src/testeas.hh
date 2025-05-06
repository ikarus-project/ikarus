// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testcommon.hh"

#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolverfactory.hh>

template <typename FE>
requires(FE::hasEAS() and FE::strainType == Ikarus::StrainTags::linear)
struct ElementTest<FE>
{
  [[nodiscard]] static auto test() {
    auto easFunctor = [](auto& f, auto& fe, auto& req, auto& affordance) {
      const auto& localView = fe.localView();
      const auto& element   = localView.element();
      constexpr int gridDim = std::remove_cvref_t<decltype(element)>::dimension;

      Dune::TestSuite t("EAS specific test");

      auto subOp = derivative(f);
      std::array<int, (gridDim == 2) ? 4 : 3> easParameters;
      if constexpr (gridDim == 2)
        easParameters = {0, 4, 5, 7};
      else
        easParameters = {0, 9, 21};

      Eigen::VectorXd oldEigenValues, newEigenValues;

      for (auto& numberOfEASParameter : easParameters) {
        fe.setEASType(numberOfEASParameter);
        auto messageIfFailed = "The number of EAS parameters are " + std::to_string(numberOfEASParameter) + ".";
        if (numberOfEASParameter == 0) {
          t.subTest(checkGradientOfElement(f, req, messageIfFailed));
          t.subTest(checkHessianOfElement(f, req, messageIfFailed));
        }
        t.subTest(checkJacobianOfElement(subOp, req, messageIfFailed));

        t.subTest(checkFEByAutoDiff(f, fe, req, affordance, messageIfFailed));

        decltype(auto) stiffnessMatrix = derivative(subOp)(req);

        auto es = Ikarus::makeIdentitySymEigenSolver<Ikarus::EigenValueSolverType::Spectra>(stiffnessMatrix);
        es.compute();
        newEigenValues = es.eigenvalues();

        t.check((newEigenValues.array() < 1e-6 * newEigenValues.norm()).count() == 3 * gridDim - 3,
                "Number of eigenvalues")
            << "We always should have " << 3 * gridDim - 3 << " zero eigenvalues, for " << 3 * gridDim - 3
            << " rigid body motions"
            << " but we counted " << (newEigenValues.array() < 1e-6 * newEigenValues.norm()).count()
            << "\nEigenValues: \n"
            << newEigenValues.transpose() << "The tolerance is " << 1e-6 * newEigenValues.norm()
            << "The number of EAS parameters is" << numberOfEASParameter << std::endl;
        if (numberOfEASParameter > 0 and numberOfEASParameter != 5)         // Q1E4 and Q1E5 are the same
          t.check((newEigenValues.array() <= oldEigenValues.array()).sum()) // Q1E4 and Q1E7 are the same in the case
                                                                            // of undistorted element
              << "More EAS parameter mean that the stiffness gets reduced. EAS parameter: " << numberOfEASParameter
              << "\noldEigenValues: \n"
              << oldEigenValues.transpose() << "\nnewEigenValues: \n"
              << newEigenValues.transpose() << std::endl;
        oldEigenValues = newEigenValues;

        const int order = 2 * (localView.tree().child(0).finiteElement().localBasis().order());
        auto rule       = Dune::QuadratureRules<double, gridDim>::rule(localView.element().type(), order);

        if (numberOfEASParameter > 0) {
          auto easVariantCopy    = fe.easVariant(); // This only test if the variant has a copy assignment operator
          const auto& easVariant = fe.easVariant();
          auto testM             = [&]<typename EAS>(const EAS& easFunction) {
            typename EAS::AnsatzType MIntegrated;
            MIntegrated.setZero();
            for (const auto& gp : rule) {
              const auto M = easFunction(gp.position());

              const double detJ = element.geometry().integrationElement(gp.position());
              MIntegrated += M * detJ * gp.weight();
            }
            t.check(MIntegrated.isZero()) << "Orthogonality condition check: The M matrix of the EAS method should be "
                                                         "zero, integrated over the domain.";
          };
          easVariant(testM);

          auto requirements = typename std::remove_cvref_t<decltype(fe)>::Requirement();
          t.checkThrow([&]() {
            calculateScalar(fe, requirements, Ikarus::ScalarAffordance::mechanicalPotentialEnergy);
          }) << "fe.calculateScalar should have failed for numberOfEASParameter > 0";
        }

        t.check(numberOfEASParameter == fe.numberOfInternalVariables())
            << "Number of EAS Parameters should be " << numberOfEASParameter << "but is "
            << fe.numberOfInternalVariables();

        t.checkThrow([&]() { fe.setEASType(100); }) << "fe.setEASType(100) should have failed";
      }
      return t;
    };
    return easFunctor;
  }
};

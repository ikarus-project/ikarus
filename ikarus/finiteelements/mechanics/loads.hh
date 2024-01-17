// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file loads.hh
 * @brief Definition of the Loads class for application of loads for finite element analysis.
 * @ingroup mechanics
 */

#pragma once

namespace Ikarus {

  /**
   * @brief Loads class represents different types of loads that can be applied.
   *
   * @ingroup mechanics
   *
   * @tparam DisplacementBasedElement The type of the displacement-based finite element.
   */
  template <typename DisplacementBasedElement>
  class Loads {
  public:
    using Traits                  = typename DisplacementBasedElement::Traits;
    using FERequirementType       = typename DisplacementBasedElement::FERequirementType;
    static constexpr int myDim    = Traits::mydim;
    static constexpr int worldDim = Traits::worlddim;

    /**
     * @brief Constructor for the Loads class.
     *
     * @param p_ele The element on which the load is applied.
     */
    explicit Loads(const DisplacementBasedElement p_ele) : ele{p_ele} {
      order = 2 * ele.localView().tree().child(0).finiteElement().localBasis().order();
    }

    /**
     * @brief Applies the distributed volume loads on the element domain.
     *
     * @tparam ScalarType The scalar type for the force vector.
     * @param par The FERequirementType object.
     * @param force The force vector to be updated.
     * @param dx Optional displacement vector.
     */
    template <typename ScalarType>
    void volume(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                const std::optional<const Eigen::VectorX<ScalarType>> dx = std::nullopt) {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = ele.displacementFunction(par, dx);
      const auto geo       = ele.localView().element().geometry();
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const Eigen::Vector<double, worldDim> fext = ele.volumeLoad(geo.global(gp.position()), lambda);
        const double intElement                    = geo.integrationElement(gp.position()) * gp.weight();
        for (size_t i = 0; i < ele.numberOfNodes; ++i) {
          const auto udCi = uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
          force.template segment<worldDim>(worldDim * i) -= udCi * fext * intElement;
        }
      }
    }

    /**
     * @brief Applies the distributed traction loads on the element boundary.
     *
     * @tparam ScalarType The scalar type for the force vector.
     * @tparam BoundaryType The type for the traction boundary.
     * @param par The FERequirementType object.
     * @param tractionBoundary The element boundary where traction load is applied.
     * @param force The force vector to be updated.
     * @param dx Optional displacement vector.
     */
    template <typename ScalarType, typename BoundaryType>
    void traction(const FERequirementType& par, const BoundaryType& tractionBoundary,
                  typename Traits::template VectorType<ScalarType> force,
                  const std::optional<const Eigen::VectorX<ScalarType>> dx = std::nullopt) {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto uFunction = ele.displacementFunction(par, dx);
      auto element         = ele.localView().element();
      for (auto&& intersection : intersections(tractionBoundary->gridView(), element)) {
        if (not tractionBoundary->contains(intersection)) continue;

        /// Integration rule along the boundary
        const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());
          const double intElement = intersection.geometry().integrationElement(curQuad.position());

          /// The value of the local function wrt the i-th coeff
          for (size_t i = 0; i < ele.numberOfNodes; ++i) {
            const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

            /// Value of the Neumann data at the current position
            auto neumannValue = ele.neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda);
            force.template segment<worldDim>(worldDim * i) -= udCi * neumannValue * curQuad.weight() * intElement;
          }
        }
      }
    }

  private:
    DisplacementBasedElement ele;
    int order{};
  };
}  // namespace Ikarus

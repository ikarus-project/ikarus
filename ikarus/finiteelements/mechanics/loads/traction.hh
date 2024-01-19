// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/fufem/boundarypatch.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/defaultfunctions.hh>

namespace Ikarus {

  /**
   * @brief Traction class represents distributed traction load that can be applied.
   *
   * @ingroup mechanics
   *
   * @tparam DisplacementBasedElement The type of the displacement-based finite element.
   * @tparam Traits Type of traits for handling finite elements.
   */
  template <typename DisplacementBasedElement, typename Traits>
  class Traction {
  public:
    using FERequirementType       = typename Traits::FERequirementType;
    using LocalView               = typename Traits::LocalView;
    using GridView                = typename Traits::GridView;
    static constexpr int myDim    = Traits::mydim;
    static constexpr int worldDim = Traits::worlddim;

    /**
     * @brief Constructor for the Loads class.
     *
     * @tparam NeumannBoundaryLoad The type for the Neumann boundary load function.
     * @param p_neumannBoundary Neumann boundary patch.
     * @param p_neumannBoundaryLoad Neumann boundary load function.
     */
    template <typename NeumannBoundaryLoad>
    explicit Traction(const BoundaryPatch<GridView>* p_neumannBoundary, NeumannBoundaryLoad p_neumannBoundaryLoad)
        : neumannBoundary{p_neumannBoundary} {
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, utils::LoadDefault>)
        neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

    /**
     * @brief Const accessor to the underlying displacement-based finite element (CRTP).
     *
     * @return Const reference to the underlying displacement-based finite element.
     */
    constexpr const DisplacementBasedElement& dbElement() const {
      return static_cast<DisplacementBasedElement const&>(*this);
    }

    /**
     * @brief Calculates the scalar energy for the given traction load.
     *
     * @tparam ScalarType The scalar type for the energy.
     *
     * @param par The FERequirementType object.
     * @param dx Optional displacement vector.
     * @return The scalar energy.
     */
    template <typename ScalarType>
    auto calculateScalar(const FERequirementType& par,
                         const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const -> ScalarType {
      const auto uFunction = dbElement().displacementFunction(par, dx);
      ScalarType energy    = 0.0;
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      auto element         = dbElement().localView().element();
      if (not neumannBoundary and not neumannBoundaryLoad) return energy;

      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), dbElement().order());

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto uVal = uFunction.evaluate(quadPos);

          // Value of the Neumann data at the current position
          auto neumannValue = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda);

          energy -= neumannValue.dot(uVal) * curQuad.weight() * integrationElement;
        }
      }
      return energy;
    }

    /**
     * @brief Calculates the force vector for the given traction load.
     *
     * @tparam ScalarType The scalar type for the force vector.
     *
     * @param par The FERequirementType object.
     * @param force The force vector to be updated.
     * @param dx Optional displacement vector.
     */
    template <typename ScalarType>
    void calculateVector(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                         const std::optional<const Eigen::VectorX<ScalarType>> dx = std::nullopt) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = dbElement().displacementFunction(par, dx);
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      auto element         = dbElement().localView().element();
      if (not neumannBoundary and not neumannBoundaryLoad) return;

      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary->contains(intersection)) continue;

        /// Integration rule along the boundary
        const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), dbElement().order());

        for (const auto& curQuad : quadLine) {
          const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());
          const double intElement = intersection.geometry().integrationElement(curQuad.position());

          /// The value of the local function wrt the i-th coeff
          for (size_t i = 0; i < dbElement().numberOfNodes(); ++i) {
            const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

            /// Value of the Neumann data at the current position
            auto neumannValue = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda);
            force.template segment<worldDim>(worldDim * i) -= udCi * neumannValue * curQuad.weight() * intElement;
          }
        }
      }
    }

    /**
     * @brief Calculates the matrix stiffness for the given traction load.
     *
     * @tparam ScalarType The scalar type for the stiffness matrix.
     *
     * @param par The FERequirementType object.
     * @param K Matrix to store the calculated stiffness.
     * @param dx Optional displacement vector.
     */
    template <typename ScalarType>
    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> K,
                         const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {}

  private:
    std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)>
        neumannBoundaryLoad;
    const BoundaryPatch<GridView>* neumannBoundary;
  };
}  // namespace Ikarus

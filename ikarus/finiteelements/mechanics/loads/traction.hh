// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/fufem/boundarypatch.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/defaultfunctions.hh>

namespace Ikarus {

/**
 * \brief Traction class represents distributed traction load that can be applied.
 *
 * \ingroup mechanics
 *
 * \tparam DFE The type of the displacement-based finite element.
 * \tparam Traits Type of traits for handling finite elements.
 */
template <typename DFE, typename Traits>
class Traction
{
public:
  using FERequirementType       = typename Traits::FERequirementType;
  using LocalView               = typename Traits::LocalView;
  using GridView                = typename Traits::GridView;
  static constexpr int myDim    = Traits::mydim;
  static constexpr int worldDim = Traits::worlddim;

  /**
   * \brief Constructor for the Loads class.
   *
   * \tparam NBL The type for the Neumann boundary load function.
   * \param neumannBoundary Neumann boundary patch.
   * \param neumannBoundaryLoad Neumann boundary load function.
   */
  template <typename NBL>
  explicit Traction(const BoundaryPatch<GridView>* neumannBoundary, NBL neumannBoundaryLoad)
      : neumannBoundary_{neumannBoundary} {
    if constexpr (!std::is_same_v<NBL, utils::LoadDefault>)
      neumannBoundaryLoad_ = neumannBoundaryLoad;

    assert(((not neumannBoundary and not neumannBoundaryLoad_) or (neumannBoundary and neumannBoundaryLoad_)) &&
           "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
  }

  /**
   * \brief Calculate the scalar value.
   *
   * Calculates the scalar value based on the given FERequirements.
   *
   * \param req The FERequirements.
   * \return The calculated scalar value.
   */
  double calculateScalar(const FERequirementType& req) const { return calculateScalarImpl<double>(req); }

  /**
   * \brief Calculate the vector associated with the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param req The FERequirementType object specifying the requirements for the calculation.
   * \param force The vector to store the calculated result.
   */
  void calculateVector(const FERequirementType& req, typename Traits::template VectorType<> force) const {
    calculateVectorImpl<double>(req, force);
  }
  /**
   * \brief Calculate the matrix associated with the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param req The FERequirementType object specifying the requirements for the calculation.
   * \param K The matrix to store the calculated result.
   */
  void calculateMatrix(const FERequirementType& req, typename Traits::template MatrixType<> K) const {
    calculateMatrixImpl<double>(req, K);
  }

protected:
  template <typename ST>
  auto calculateScalarImpl(const FERequirementType& par,
                           const std::optional<const Eigen::VectorX<ST>>& dx = std::nullopt) const -> ST {
    if (not neumannBoundary_ and not neumannBoundaryLoad_)
      return 0.0;
    ST energy            = 0.0;
    const auto uFunction = dbElement().displacementFunction(par, dx);
    const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
    auto& element        = dbElement().localView().element();

    for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
      if (not neumannBoundary_->contains(intersection))
        continue;

      const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), dbElement().order());

      for (const auto& curQuad : quadLine) {
        // Local position of the quadrature point
        const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());

        const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

        // The value of the local function
        const auto uVal = uFunction.evaluate(quadPos);

        // Value of the Neumann data at the current position
        auto neumannValue = neumannBoundaryLoad_(intersection.geometry().global(curQuad.position()), lambda);

        energy -= neumannValue.dot(uVal) * curQuad.weight() * integrationElement;
      }
    }
    return energy;
  }

  template <typename ST>
  void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ST> force,
                           const std::optional<const Eigen::VectorX<ST>> dx = std::nullopt) const {
    if (not neumannBoundary_ and not neumannBoundaryLoad_)
      return;
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = dbElement().displacementFunction(par, dx);
    const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
    auto& element        = dbElement().localView().element();

    for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
      if (not neumannBoundary_->contains(intersection))
        continue;

      /// Integration rule along the boundary
      const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), dbElement().order());

      for (const auto& curQuad : quadLine) {
        const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());
        const double intElement = intersection.geometry().integrationElement(curQuad.position());

        /// The value of the local function wrt the i-th coeff
        for (size_t i = 0; i < dbElement().numberOfNodes(); ++i) {
          const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

          /// Value of the Neumann data at the current position
          auto neumannValue = neumannBoundaryLoad_(intersection.geometry().global(curQuad.position()), lambda);
          force.template segment<worldDim>(worldDim * i) -= udCi * neumannValue * curQuad.weight() * intElement;
        }
      }
    }
  }

  template <typename ST>
  void calculateMatrixImpl(const FERequirementType& par, typename Traits::template MatrixType<> K,
                           const std::optional<const Eigen::VectorX<ST>>& dx = std::nullopt) const {}

private:
  std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)>
      neumannBoundaryLoad_;
  const BoundaryPatch<GridView>* neumannBoundary_;

  /**
   * \brief Const accessor to the underlying displacement-based finite element (CRTP).
   *
   * \return Const reference to the underlying displacement-based finite element.
   */
  constexpr const DFE& dbElement() const { return static_cast<const DFE&>(*this); }
};
} // namespace Ikarus

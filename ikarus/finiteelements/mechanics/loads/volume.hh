// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/defaultfunctions.hh>

namespace Ikarus {

/**
 * \brief Volume class represents distributed volume load that can be applied.
 *
 * \ingroup mechanics
 *
 * \tparam DFE The type of the displacement-based finite element.
 * \tparam Traits Type of traits for handling finite elements.
 */
template <typename DFE, typename Traits>
class Volume
{
public:
  using FERequirementType       = typename Traits::FERequirementType;
  using LocalView               = typename Traits::LocalView;
  static constexpr int worldDim = Traits::worlddim;

  /**
   * \brief Constructor for the Loads class.
   *
   * \tparam VL The type for the volume load function.
   *
   * \param volumeLoad Volume load function.
   */
  template <typename VL>
  explicit Volume(VL volumeLoad = {}) {
    if constexpr (!std::is_same_v<VL, utils::LoadDefault>)
      volumeLoad_ = volumeLoad;
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
    if (not volumeLoad_)
      return 0.0;
    ST energy            = 0.0;
    const auto uFunction = dbElement().displacementFunction(par, dx);
    const auto geo       = dbElement().geometry();
    const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);

    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const auto uVal                      = uFunction.evaluate(gpIndex);
      Eigen::Vector<double, worldDim> fext = volumeLoad_(geo.global(gp.position()), lambda);
      energy -= uVal.dot(fext) * geo.integrationElement(gp.position()) * gp.weight();
    }
    return energy;
  }

  template <typename ST>
  void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ST> force,
                           const std::optional<const Eigen::VectorX<ST>>& dx = std::nullopt) const {
    if (not volumeLoad_)
      return;
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = dbElement().displacementFunction(par, dx);
    const auto geo       = dbElement().geometry();
    const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);

    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const Eigen::Vector<double, worldDim> fext = volumeLoad_(geo.global(gp.position()), lambda);
      const double intElement                    = geo.integrationElement(gp.position()) * gp.weight();
      for (size_t i = 0; i < dbElement().numberOfNodes(); ++i) {
        const auto udCi = uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
        force.template segment<worldDim>(worldDim * i) -= udCi * fext * intElement;
      }
    }
  }

  template <typename ST>
  void calculateMatrixImpl(const FERequirementType& par, typename Traits::template MatrixType<> K,
                           const std::optional<const Eigen::VectorX<ST>>& dx = std::nullopt) const {}

private:
  std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)> volumeLoad_;

  /**
   * \brief Const accessor to the underlying displacement-based finite element (CRTP).
   *
   * \return Const reference to the underlying displacement-based finite element.
   */
  constexpr const DFE& dbElement() const { return static_cast<const DFE&>(*this); }
};
} // namespace Ikarus

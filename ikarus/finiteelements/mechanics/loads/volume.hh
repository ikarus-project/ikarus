// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/defaultfunctions.hh>

namespace Ikarus {

  /**
   * @brief Volume class represents distributed volume load that can be applied.
   *
   * @ingroup mechanics
   *
   * @tparam DisplacementBasedElement The type of the displacement-based finite element.
   * @tparam Traits Type of traits for handling finite elements.
   */
  template <typename DisplacementBasedElement, typename Traits>
  class Volume {
  public:
    using FERequirementType       = typename Traits::FERequirementType;
    using LocalView               = typename Traits::LocalView;
    static constexpr int worldDim = Traits::worlddim;

    /**
     * @brief Constructor for the Loads class.
     *
     * @tparam VolumeLoad The type for the volume load function.
     *
     * @param p_volumeLoad Volume load function.
     */
    template <typename VolumeLoad>
    explicit Volume(VolumeLoad p_volumeLoad = {}) {
      if constexpr (!std::is_same_v<VolumeLoad, utils::LoadDefault>) volumeLoad = p_volumeLoad;
    }

    /**
     * @brief Calculate the scalar value.
     *
     * Calculates the scalar value based on the given FERequirements.
     *
     * @param req The FERequirements.
     * @return The calculated scalar value.
     */
    double calculateScalar(const FERequirementType& req) const { return calculateScalarImpl<double>(req); }

    /**
     * @brief Calculate the vector associated with the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the calculation.
     * @param req The FERequirementType object specifying the requirements for the calculation.
     * @param force The vector to store the calculated result.
     */
    void calculateVector(const FERequirementType& req, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(req, force);
    }
    /**
     * @brief Calculate the matrix associated with the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the calculation.
     * @param req The FERequirementType object specifying the requirements for the calculation.
     * @param K The matrix to store the calculated result.
     */
    void calculateMatrix(const FERequirementType& req, typename Traits::template MatrixType<> K) const {
      calculateMatrixImpl<double>(req, K);
    }

  protected:
    template <typename ScalarType>
    auto calculateScalarImpl(const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                                           = std::nullopt) const -> ScalarType {
      if (not volumeLoad) return 0.0;
      ScalarType energy    = 0.0;
      const auto uFunction = dbElement().displacementFunction(par, dx);
      const auto geo       = dbElement().geometry();
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto uVal                      = uFunction.evaluate(gpIndex);
        Eigen::Vector<double, worldDim> fext = volumeLoad(geo.global(gp.position()), lambda);
        energy -= uVal.dot(fext) * geo.integrationElement(gp.position()) * gp.weight();
      }
      return energy;
    }

    template <typename ScalarType>
    void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      if (not volumeLoad) return;
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = dbElement().displacementFunction(par, dx);
      const auto geo       = dbElement().geometry();
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const Eigen::Vector<double, worldDim> fext = volumeLoad(geo.global(gp.position()), lambda);
        const double intElement                    = geo.integrationElement(gp.position()) * gp.weight();
        for (size_t i = 0; i < dbElement().numberOfNodes(); ++i) {
          const auto udCi = uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
          force.template segment<worldDim>(worldDim * i) -= udCi * fext * intElement;
        }
      }
    }

    template <typename ScalarType>
    void calculateMatrixImpl(const FERequirementType& par, typename Traits::template MatrixType<> K,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {}

  private:
    std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)>
        volumeLoad;

    /**
     * @brief Const accessor to the underlying displacement-based finite element (CRTP).
     *
     * @return Const reference to the underlying displacement-based finite element.
     */
    constexpr const DisplacementBasedElement& dbElement() const {
      return static_cast<DisplacementBasedElement const&>(*this);
    }
  };
}  // namespace Ikarus

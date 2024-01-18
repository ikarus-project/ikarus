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
     * @param p_dispFE Displacement-based finite element.
     * @param p_volumeLoad Volume load function.
     */
    template <typename VolumeLoad>
    explicit Volume(const DisplacementBasedElement& p_dispFE, VolumeLoad p_volumeLoad = {})
        : dispFE{std::make_shared<DisplacementBasedElement>(p_dispFE)} {
      if constexpr (!std::is_same_v<VolumeLoad, utils::LoadDefault>) volumeLoad = p_volumeLoad;
    }

    /**
     * @brief Calculates the scalar energy for the given FERequirementType.
     *
     * @param par The FERequirementType object.
     * @return The scalar energy.
     */
    template <typename ScalarType>
    auto calculateScalar(const FERequirementType& par,
                         const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const -> ScalarType {
      const auto uFunction = dispFE->displacementFunction(par, dx);
      const auto geo       = dispFE->geometry();
      ScalarType energy    = 0.0;
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      if (not(volumeLoad)) return energy;

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto uVal                      = uFunction.evaluate(gpIndex);
        Eigen::Vector<double, worldDim> fext = volumeLoad(geo->global(gp.position()), lambda);
        energy -= uVal.dot(fext) * geo->integrationElement(gp.position()) * gp.weight();
      }
      return energy;
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
    void calculateVector(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                         const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto& localView = dispFE->localView();
      std::cout << "Indices in Volume class\n";
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 4; ++j) {
          const auto localIndex = localView.tree().child(i).localIndex(j);
          std::cout << "localIndex\t" << localIndex << std::endl;
          const auto globalIndex = localView.index(localIndex)[0];
          std::cout << "globalIndex\t" << globalIndex << std::endl;
        }
      std::cout << "localView().size()\t" << localView.size() << std::endl;
      std::cout << "finiteElement.size()\t" << localView.tree().child(0).finiteElement().size() << std::endl;
      const auto uFunction = dispFE->displacementFunction(par, dx);
      const auto geo       = dispFE->geometry();
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      if (not(volumeLoad)) return;

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const Eigen::Vector<double, worldDim> fext = volumeLoad(geo->global(gp.position()), lambda);
        const double intElement                    = geo->integrationElement(gp.position()) * gp.weight();
        for (size_t i = 0; i < dispFE->numberOfNodes(); ++i) {
          const auto udCi = uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
          force.template segment<worldDim>(worldDim * i) -= udCi * fext * intElement;
        }
      }
    }

    /**
     * @brief Calculates the matrix stiffness for the given FERequirementType.
     *
     * @param par The FERequirementType object.
     * @param K Matrix to store the calculated stiffness.
     */
    template <typename ScalarType>
    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> K,
                         const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {}

  private:
    const std::shared_ptr<DisplacementBasedElement> dispFE;
    std::function<Eigen::Vector<double, worldDim>(const Dune::FieldVector<double, worldDim>&, const double&)>
        volumeLoad;
  };
}  // namespace Ikarus

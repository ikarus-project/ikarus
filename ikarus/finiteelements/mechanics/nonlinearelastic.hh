// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file nonlinearelastic.hh
 * @brief Definition of the NonLinearElastic class for finite element mechanics computations.
 * @ingroup  mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS
#  include <dune/fufem/boundarypatch.hh>
#  include <dune/geometry/quadraturerules.hh>
#  include <dune/geometry/type.hh>
#  include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#  include <dune/localfefunctions/expressions/greenLagrangeStrains.hh>
#  include <dune/localfefunctions/impl/standardLocalFunction.hh>
#  include <dune/localfefunctions/manifolds/realTuple.hh>

#  include <ikarus/finiteelements/febases/powerbasisfe.hh>
#  include <ikarus/finiteelements/ferequirements.hh>
#  include <ikarus/finiteelements/fetraits.hh>
#  include <ikarus/finiteelements/mechanics/materials/tags.hh>
#  include <ikarus/finiteelements/physicshelper.hh>
#  include <ikarus/utils/defaultfunctions.hh>
#  include <ikarus/utils/eigendunetransformations.hh>
#  include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

  /**
   * @brief NonLinearElastic class represents a non-linear elastic finite element.
   *
   * @ingroup mechanics
   *
   * @tparam Basis_ The basis type for the finite element.
   * @tparam Material_ The material type for the finite element.
   * @tparam FERequirements_ The requirements for the finite element.
   * @tparam useEigenRef A boolean flag indicating whether to use Eigen references.
   */
  template <typename Basis_, typename Material_, typename FERequirements_ = FERequirements<>, bool useEigenRef = false>
  class NonLinearElastic : public PowerBasisFE<typename Basis_::FlatBasis> {
  public:
    using Basis                      = Basis_;
    using Material                   = Material_;
    using FlatBasis                  = typename Basis::FlatBasis;
    using BasePowerFE                = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType          = FERequirements_;
    using ResultRequirementsType     = ResultRequirements<FERequirementType>;
    using LocalView                  = typename FlatBasis::LocalView;
    using Element                    = typename LocalView::Element;
    using Geometry                   = typename Element::Geometry;
    using GridView                   = typename FlatBasis::GridView;
    using Traits                     = TraitsFromLocalView<LocalView, useEigenRef>;
    static constexpr int myDim       = Traits::mydim;
    static constexpr auto strainType = StrainTags::greenLagrangian;

    /**
     * @brief Constructor for the NonLinearElastic class.
     *
     * @tparam VolumeLoad The type for the volume load function.
     * @tparam NeumannBoundaryLoad The type for the Neumann boundary load function.
     * @param globalBasis The global basis for the finite element.
     * @param element The element for which the finite element is constructed.
     * @param p_mat The material for the non-linear elastic element.
     * @param p_volumeLoad Volume load function (default is LoadDefault).
     * @param p_neumannBoundary Neumann boundary patch (default is nullptr).
     * @param p_neumannBoundaryLoad Neumann boundary load function (default is LoadDefault).
     */
    template <typename VolumeLoad = utils::LoadDefault, typename NeumannBoundaryLoad = utils::LoadDefault>
    NonLinearElastic(const Basis& globalBasis, const typename LocalView::Element& element, const Material& p_mat,
                     VolumeLoad p_volumeLoad = {}, const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
                     NeumannBoundaryLoad p_neumannBoundaryLoad = {})
        : BasePowerFE(globalBasis.flat(), element), neumannBoundary{p_neumannBoundary}, mat{p_mat} {
      this->localView().bind(element);
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes     = fe.size();
      dispAtNodes.resize(numberOfNodes);
      order      = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      localBasis = Dune::CachedLocalBasis(this->localView().tree().child(0).finiteElement().localBasis());
      if constexpr (requires { this->localView().element().impl().getQuadratureRule(order); })
        if (this->localView().element().impl().isTrimmed())
          localBasis.bind(this->localView().element().impl().getQuadratureRule(order), Dune::bindDerivatives(0, 1));
        else
          localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                          Dune::bindDerivatives(0, 1));
      else
        localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                        Dune::bindDerivatives(0, 1));

      if constexpr (!std::is_same_v<VolumeLoad, utils::LoadDefault>) volumeLoad = p_volumeLoad;
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, utils::LoadDefault>)
        neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

    /**
     * @brief Get the displacement function for the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the displacement function.
     * @param par The FERequirementType object.
     * @param dx Optional displacement vector.
     * @return A StandardLocalFunction representing the displacement function.
     */
    template <typename ScalarType = double>
    auto displacementFunction(const FERequirementType& par,
                              const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);

      Dune::BlockVector<Dune::RealTuple<ScalarType, Traits::dimension>> disp(dispAtNodes.size());
      // If optional displacement vector is provided, apply it
      if (dx)
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < myDim; ++k2)
            disp[i][k2] = dx.value()[i * myDim + k2]
                          + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      else
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < myDim; ++k2)
            disp[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];

      auto geo = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          this->localView().element().geometry());
      Dune::StandardLocalFunction uFunction(localBasis, disp, geo);
      return uFunction;
    }

    /**
     * @brief The strain function for the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the strain function.
     * @param par The FERequirementType object.
     * @param dx Optional displacement vector.
     * @return The strain function calculated using greenLagrangeStrains.
     */
    template <typename ScalarType = double>
    auto strainFunction(const FERequirementType& par,
                        const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      return greenLagrangeStrains(displacementFunction(par, dx));
    }

    /**
     * @brief Get the material tangent for the given strain.
     *
     * @tparam ScalarType The scalar type for the material and strain.
     * @tparam strainDim The dimension of the strain vector.
     * @tparam voigt Flag indicating whether to use Voigt notation.
     * @param strain The strain vector.
     * @return The material tangent calculated using the material's tangentModuli function.
     */
    template <typename ScalarType, int strainDim, bool voigt = true>
    auto getMaterialTangent(const Eigen::Vector<ScalarType, strainDim>& strain) const {
      if constexpr (std::is_same_v<ScalarType, double>)
        return mat.template tangentModuli<strainType, voigt>(strain);
      else {
        decltype(auto) matAD = mat.template rebind<ScalarType>();
        return matAD.template tangentModuli<strainType, voigt>(strain);
      }
    }

    /**
     * @brief Get the internal energy for the given strain.
     *
     * @tparam ScalarType The scalar type for the material and strain.
     * @tparam strainDim The dimension of the strain vector.
     * @param strain The strain vector.
     * @return The internal energy calculated using the material's storedEnergy function.
     */
    template <typename ScalarType, int strainDim>
    auto getInternalEnergy(const Eigen::Vector<ScalarType, strainDim>& strain) const {
      if constexpr (std::is_same_v<ScalarType, double>)
        return mat.template storedEnergy<strainType>(strain);
      else {
        decltype(auto) matAD = mat.template rebind<ScalarType>();
        return matAD.template storedEnergy<strainType>(strain);
      }
    }

    /**
     * @brief Get the stress for the given strain.
     *
     * @tparam ScalarType The scalar type for the material and strain.
     * @tparam strainDim The dimension of the strain vector.
     * @tparam voigt A boolean indicating whether to use the Voigt notation for stress.
     * @param strain The strain vector.
     * @return The stress vector calculated using the material's stresses function.
     */
    template <typename ScalarType, int strainDim, bool voigt = true>
    auto getStress(const Eigen::Vector<ScalarType, strainDim>& strain) const {
      if constexpr (std::is_same_v<ScalarType, double>)
        return mat.template stresses<strainType, voigt>(strain);
      else {
        decltype(auto) matAD = mat.template rebind<ScalarType>();
        return matAD.template stresses<strainType, voigt>(strain);
      }
    }

    /**
     * @brief Calculate the scalar value associated with the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the calculation.
     * @param par The FERequirementType object specifying the requirements for the calculation.
     * @return The calculated scalar value.
     */
    double calculateScalar(const FERequirementType& par) const { return calculateScalarImpl<double>(par); }

    /**
     * @brief Calculate the vector associated with the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the calculation.
     * @param par The FERequirementType object specifying the requirements for the calculation.
     * @param force The vector to store the calculated result.
     */
    void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(par, force);
    }

    /**
     * @brief Calculate the matrix associated with the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the calculation.
     * @param par The FERequirementType object specifying the requirements for the calculation.
     * @param K The matrix to store the calculated result.
     */
    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> K) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto eps = strainFunction(par);
      const auto geo = this->localView().element().geometry();

      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        const auto EVoigt       = (eps.evaluate(gpIndex, on(gridElement))).eval();
        const auto C            = getMaterialTangent(EVoigt);
        const auto stresses     = getStress(EVoigt);
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          for (size_t j = 0; j < numberOfNodes; ++j) {
            const auto bopJ = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), on(gridElement));
            const auto kgIJ = eps.evaluateDerivative(gpIndex, wrt(coeff(i, j)), along(stresses), on(gridElement));
            K.template block<myDim, myDim>(i * myDim, j * myDim) += (bopI.transpose() * C * bopJ + kgIJ) * intElement;
          }
        }
      }
    }

    /**
     * @brief Calculate specified results at a given local position.
     *
     * @param req The ResultRequirementsType object specifying the required results.
     * @param local The local position for which results are to be calculated.
     * @param result The ResultTypeMap object to store the calculated results.
     */
    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto uFunction = displacementFunction(req.getFERequirements());
      const auto H         = uFunction.evaluateDerivative(local, Dune::wrt(spatialAll), Dune::on(gridElement));
      const auto E         = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
      const auto EVoigt    = toVoigt(E);
      auto PK2             = mat.template stresses<StrainTags::greenLagrangian>(EVoigt);

      if (req.isResultRequested(ResultType::PK2Stress))
        result.insertOrAssignResult(ResultType::PK2Stress, PK2);
      else
        DUNE_THROW(Dune::NotImplemented, "The requested result type is NOT implemented.");
    }

    Dune::CachedLocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Dune::FieldVector<double, Traits::worlddim>&,
                                                          const double&)>
        volumeLoad;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Dune::FieldVector<double, Traits::worlddim>&,
                                                          const double&)>
        neumannBoundaryLoad;
    const BoundaryPatch<GridView>* neumannBoundary;
    Material mat;
    mutable Dune::BlockVector<Dune::RealTuple<double, Traits::dimension>> dispAtNodes;
    size_t numberOfNodes{0};
    int order{};

  protected:
    template <typename ScalarType>
    auto calculateScalarImpl(const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                                           = std::nullopt) const -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = displacementFunction(par, dx);
      const auto eps       = strainFunction(par, dx);
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto geo       = this->localView().element().geometry();
      ScalarType energy    = 0.0;

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto EVoigt         = (eps.evaluate(gpIndex, on(gridElement))).eval();
        const auto internalEnergy = getInternalEnergy(EVoigt);
        energy += internalEnergy * geo.integrationElement(gp.position()) * gp.weight();
      }

      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
          const auto u                                       = uFunction.evaluate(gpIndex);
          const Eigen::Vector<double, Traits::worlddim> fExt = volumeLoad(geo.global(gp.position()), lambda);
          energy -= u.dot(fExt) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e., neumann boundary
      if (not neumannBoundary and not neumannBoundaryLoad) return energy;

      const auto& element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary or not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, Traits::mydim - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, Traits::mydim>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double intElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto u = uFunction.evaluate(quadPos);

          // Value of the Neumann data at the current position
          const auto neumannValue = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda);
          energy -= neumannValue.dot(u) * curQuad.weight() * intElement;
        }
      }
      return energy;
    }

    template <typename ScalarType>
    void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto eps     = strainFunction(par, dx);
      const auto geo     = this->localView().element().geometry();

      // Internal forces
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        const auto EVoigt       = (eps.evaluate(gpIndex, on(gridElement))).eval();
        const auto stresses     = getStress(EVoigt);
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          force.template segment<myDim>(myDim * i) += bopI.transpose() * stresses * intElement;
        }
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        const auto u = displacementFunction(par, dx);
        for (const auto& [gpIndex, gp] : u.viewOverIntegrationPoints()) {
          const double intElement                            = geo.integrationElement(gp.position()) * gp.weight();
          const Eigen::Vector<double, Traits::worlddim> fExt = volumeLoad(geo.global(gp.position()), lambda);
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<myDim>(myDim * i) -= udCi * fExt * intElement;
          }
        }
      }

      // External forces, boundary forces, i.e., at the Neumann boundary
      if (not neumannBoundary and not neumannBoundaryLoad) return;

      const auto u        = displacementFunction(par, dx);
      const auto& element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary->contains(intersection)) continue;

        // Integration rule along the boundary
        const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double intElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function wrt the i-th coefficient
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(quadPos, wrt(coeff(i)));

            // Value of the Neumann data at the current position
            const auto neumannValue = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda);
            force.template segment<myDim>(myDim * i) -= udCi * neumannValue * curQuad.weight() * intElement;
          }
        }
      }
    }
  };
}  // namespace Ikarus

#else
#  error NonLinearElastic depends on dune-localfefunctions, which is not included
#endif

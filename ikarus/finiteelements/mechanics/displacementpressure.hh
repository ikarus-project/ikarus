// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file displacementpressure.hh
 * \brief Definition of the DisplacementPressure class for finite element mechanics computations.
 * \ingroup mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/fufem/boundarypatch.hh>
  #include <dune/geometry/quadraturerules.hh>
  #include <dune/geometry/type.hh>
  #include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
  #include <dune/localfefunctions/expressions/greenLagrangeStrains.hh>
  #include <dune/localfefunctions/impl/standardLocalFunction.hh>
  #include <dune/localfefunctions/manifolds/realTuple.hh>

  #include <ikarus/finiteelements/febase.hh>
  #include <ikarus/finiteelements/fehelper.hh>
  #include <ikarus/finiteelements/ferequirements.hh>
  #include <ikarus/finiteelements/mechanics/loads.hh>
  #include <ikarus/finiteelements/mechanics/materials/decomposehyperelastic.hh>
  #include <ikarus/finiteelements/mechanics/materials/hyperelastic/interface.hh>
  #include <ikarus/finiteelements/mechanics/materials/tags.hh>
  #include <ikarus/finiteelements/mechanics/strainenhancements/easfunctions.hh>
  #include <ikarus/finiteelements/physicshelper.hh>
  #include <ikarus/utils/defaultfunctions.hh>
  #include <ikarus/utils/eigendunetransformations.hh>
  #include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

namespace Impl {

  /**
   * \brief Constitutive driver for the displacement pressure element using a split between a volumetric and deviatoric
   * function.
   *
   * \tparam DEV a Ikarus::Material representing the deviatoric part of the material
   * \tparam VOL a Ikarus::Material representing the volumetric part of the material
   */
  template <Concepts::Material DEV, Concepts::Material VOL>
  struct DPConstitutiveDriver
  {
    static constexpr auto strainType = StrainTags::greenLagrangian;

    auto Uhat(const auto& Estrain) const { return vol_.template storedEnergy<strainType>(Estrain); }

    auto storedEnergy(const auto& Estrain, const auto& p) {
      auto Wdev = dev_.template storedEnergy<strainType>(Estrain);
      auto uhat = Uhat(Estrain);
      return Wdev + p * uhat - (1.0 / (2.0 * kappa_)) * Dune::power(p, 2);
    }

    auto firstDerivatives(const auto& Estrain, const auto& p) const {
      const auto Sdev = dev_.template stresses<strainType, true>(Estrain);
      const auto dUdE = vol_.template stresses<strainType, true>(Estrain);
      const auto S    = (Sdev + p * dUdE).eval();
      return std::make_tuple(dUdE, S);
    }

    auto secondDerivative(const auto& Estrain, const auto& p) const {
      const auto Cdev   = dev_.template tangentModuli<strainType, true>(Estrain);
      const auto d2UdE2 = vol_.template tangentModuli<strainType, true>(Estrain);
      const auto C      = (Cdev + p * d2UdE2).eval();
      return C;
    }

    double kappa() const { return kappa_; }

    DEV deviatoricMaterial() const { return dev_; }
    VOL volumetricMaterial() const { return vol_; }

    DPConstitutiveDriver(const DEV& dev, const VOL& vol, double kappa)
        : dev_(dev),
          vol_(vol),
          kappa_(kappa) {}

    /**
     * \brief Rebinds the material driver to a different scalar type.
     * \tparam STO The target scalar type.
     * \return ConstitutiveDriver<ScalarTypeOther> The rebound ConstitutiveDriver.
     */
    template <typename STO>
    auto rebind() const {
      auto reboundDEV = dev_.template rebind<STO>();
      auto reboundVOL = vol_.template rebind<STO>();
      return DPConstitutiveDriver<decltype(reboundDEV), decltype(reboundVOL)>(reboundDEV, reboundVOL, kappa_);
    }

  private:
    DEV dev_;
    VOL vol_;
    double kappa_;
  };

} // namespace Impl
template <typename PreFE, typename FE, typename PRE>
class DisplacementPressure;

/**
 * \brief A PreFE struct for displacement-pressure elements.
 * \tparam MAT Type of the material.
 */
template <Concepts::Material MAT>
struct DisplacementPressurePre
{
  using Material = MAT;
  MAT material;

  template <typename PreFE, typename FE>
  using Skill = DisplacementPressure<PreFE, FE, DisplacementPressurePre>;
};

/**
 * \brief DisplacementPressure class represents a displacement-pressure finite element.
 *
 * \ingroup mechanics
 *
 * \tparam PreFE The type of the total pre finite element.
 * \tparam FE The type of the finite element.
 * \tparam PRE The type of the displacement-pressure pre finite element.
 */
template <typename PreFE, typename FE, typename PRE>
class DisplacementPressure : public ResultTypeBase<ResultTypes::PK2Stress, ResultTypes::PK2StressFull>
{
public:
  using Traits      = PreFE::Traits;
  using Basis       = typename Traits::Basis;
  using FlatBasis   = typename Traits::FlatBasis;
  using Requirement = FERequirements<FESolutions::displacement, FEParameter::loadfactor>;
  using LocalView   = typename Traits::LocalView;
  using Geometry    = typename Traits::Geometry;
  using GridView    = typename Traits::GridView;
  using Element     = typename Traits::Element;
  using Pre         = PRE;

  using LocalBasisTypeU =
      decltype(std::declval<LocalView>().tree().child(Dune::Indices::_0).child(0).finiteElement().localBasis());
  using LocalBasisTypeP =
      decltype(std::declval<LocalView>().tree().child(Dune::Indices::_1).finiteElement().localBasis());

  template <typename ST>
  using VectorXOptRef = std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>;

  static constexpr int myDim       = Traits::mydim;
  static constexpr int strainDim   = myDim * (myDim + 1) / 2;
  static constexpr auto strainType = StrainTags::greenLagrangian;
  static constexpr auto stressType = StressTags::PK2;

  template <template <typename, int, int> class RT>
  using RTWrapperType = ResultWrapper<RT<typename Traits::ctype, myDim, Traits::worlddim>, ResultShape::Vector>;

  template <typename ST = double>
  using StrainType = Eigen::Vector<ST, strainDim>;

  template <typename ST = double>
  using BopType = Eigen::Matrix<ST, strainDim, myDim>;

  template <typename ST = double>
  using KgType = Eigen::Matrix<ST, myDim, myDim>;

  using MaterialType = typename Pre::Material;

  template <typename MAT>
  using DecomposedDevType = typename Materials::DecomposedMaterialTypes<MaterialType>::DevType;

  template <typename MAT>
  using DecomposedVolType = typename Materials::DecomposedMaterialTypes<MaterialType>::VolType;

  using DPConstitutiveDriverType =
      Impl::DPConstitutiveDriver<DecomposedDevType<MaterialType>, DecomposedVolType<MaterialType>>;

  /**
   * \brief Constructor for the DisplacementPressure class.
   * \param pre The pre fe
   */
  explicit DisplacementPressure(const Pre& pre)
      : constitutiveDriver_([](const auto& mat) {
          auto [dev, vol, devParams, volParams] = Materials::decomposeHyperelasticAndGetMaterialParameters(mat);
          return DPConstitutiveDriverType(dev, vol, volParams);
        }(pre.material)) {
    static_assert(Pre::Material::isHyperelastic, "DisplacementPressure is only implemented for the hyperelastic case.");
  }

protected:
  /**
   * \brief A helper function to bind the local view to the element.
   */
  void bindImpl() {
    const auto& localView = underlying().localView();
    const auto& element   = localView.element();
    auto& firstChild      = localView.tree().child(Dune::Indices::_0);
    const auto& firstFE   = firstChild.child(0).finiteElement();

    auto& secondChild    = localView.tree().child(Dune::Indices::_1);
    const auto& secondFE = secondChild.finiteElement();

    geo_           = std::make_shared<const Geometry>(element.geometry());
    numberOfNodes_ = firstFE.size();
    order_         = 2 * (firstFE.localBasis().order());
    localBasisU_   = Dune::CachedLocalBasis(firstFE.localBasis());
    localBasisP_   = Dune::CachedLocalBasis(secondFE.localBasis());

    localBasisU_.bind(Dune::QuadratureRules<double, myDim>::rule(element.type(), order_), Dune::bindDerivatives(0, 1));
    localBasisP_.bind(Dune::QuadratureRules<double, myDim>::rule(element.type(), order_), Dune::bindDerivatives(0, 1));
  }

public:
  /**
   * \brief Get the displacement function for the given Requirement.
   *
   * \tparam ScalarType The scalar type for the displacement function.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return A StandardLocalFunction representing the displacement function.
   */
  template <typename ScalarType = double>
  auto displacementFunction(const Requirement& par, const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    const auto& d = par.globalSolution();
    auto disp =
        Ikarus::FEHelper::localSolutionBlockVectorComposite<Traits>(d, underlying().localView(), Dune::Indices::_0, dx);
    Dune::StandardLocalFunction uFunction(localBasisU_, disp, geo_);
    return uFunction;
  }

  /**
   * \brief Get the pressure function for the given Requirement.
   *
   * \tparam ScalarType The scalar type for the pressure function.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return A StandardLocalFunction representing the pressure function.
   */
  template <typename ScalarType = double>
  auto pressureFunction(const Requirement& par, const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    const auto& d = par.globalSolution();
    auto pressure =
        Ikarus::FEHelper::localSolutionBlockVectorScalar<Traits>(d, underlying().localView(), Dune::Indices::_1, dx);

    Dune::StandardLocalFunction pFunction(localBasisP_, pressure, geo_);
    return pFunction;
  }

  /**
   * \brief The strain function for the given Requirement.
   *
   * \tparam ScalarType The scalar type for the strain function.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return The strain function calculated using greenLagrangeStrains.
   */
  template <typename ScalarType = double>
  inline auto strainFunction(const Requirement& par, const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    return Dune::greenLagrangeStrains(displacementFunction(par, dx));
  }

  const Geometry& geometry() const { return *geo_; }
  [[nodiscard]] size_t numberOfNodes() const { return numberOfNodes_; }
  [[nodiscard]] int order() const { return order_; }
  const Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisTypeU>>& localBasis() const { return localBasisU_; }
  const Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisTypeP>>& localBasisP() const { return localBasisP_; }

  template <typename ScalarType = double>
  decltype(auto) constitutiveDriver() const {
    if constexpr (Concepts::AutodiffScalar<ScalarType>)
      return constitutiveDriver_.template rebind<ScalarType>();
    else
      return constitutiveDriver_;
  }

public:
  /**
   * \brief Get a lambda function that evaluates the requested result type for a given strain (in Voigt notation)
   * and pressure.
   * \tparam RT The type representing the requested result.
   * \return A lambda function that evaluates the requested result type for a given strain (in Voigt notation).
   */
  template <template <typename, int, int> class RT>
  requires(supportsResultType<RT>())
  auto resultFunction() const {
    return [&](const Eigen::Vector<double, strainDim>& strainInVoigt, double p) {
      if constexpr (isSameResultType<RT, ResultTypes::PK2Stress> or isSameResultType<RT, ResultTypes::PK2StressFull>) {
        auto [dev, vol] = [&]() {
          if constexpr (isSameResultType<RT, ResultTypes::PK2StressFull> and PRE::Material::isReduced)
            return std::make_pair(constitutiveDriver_.deviatoricMaterial().underlying(),
                                  constitutiveDriver_.volumetricMaterial().underlying());
          else
            return std::make_pair(constitutiveDriver_.deviatoricMaterial(), constitutiveDriver_.volumetricMaterial());
        }();

        return RTWrapperType<RT>{
            dev.template stresses<strainType>(enlargeIfReduced<typename PRE::Material>(strainInVoigt)) +
            p * vol.template stresses<strainType>(enlargeIfReduced<typename PRE::Material>(strainInVoigt))};
      }
    };
  }

  /**
   * \brief Calculates a requested result at a specific local position.
   *
   * \param req The Requirement object holding the global solution.
   * \param local Local position vector.
   * \return calculated result
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  requires(supportsResultType<RT>())
  auto calculateAtImpl(const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<1>) const {
    using namespace Dune::DerivativeDirections;

    if constexpr (isSameResultType<RT, ResultTypes::PK2Stress> or isSameResultType<RT, ResultTypes::PK2StressFull>) {
      const auto uFunction = displacementFunction(req);
      const auto pFunction = pressureFunction(req);

      const auto rFunction = resultFunction<RT>();
      const auto H         = uFunction.evaluateDerivative(local, Dune::wrt(spatialAll), Dune::on(gridElement));
      const auto p         = pFunction.evaluate(local, Dune::on(gridElement)).eval();
      const auto E         = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();

      return rFunction(toVoigt(E), p[0]);
    }
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }
  std::shared_ptr<const Geometry> geo_;
  Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisTypeU>> localBasisU_;
  Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisTypeP>> localBasisP_;
  DPConstitutiveDriverType constitutiveDriver_;
  size_t numberOfNodes_{0};
  int order_{};

public:
  /**
   * \brief Get a lambda function that evaluates the geometric part of the stiffness matrix (Kg) for a given integration
   * point and its index.
   *
   * \tparam ST The scalar type for the material and strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \param K The matrix to store the calculated result.
   * \return A lambda function that evaluates the stiffness matrix for a given strain, integration point and its index.
   */
  template <typename ST>
  auto geometricStiffnessMatrixFunction(const Requirement& par, typename Traits::template MatrixType<ST>& K,
                                        const VectorXOptRef<ST>& dx = std::nullopt) const {
    return [&](const FE::template KgType<ST>& kgIJ, const int I, const int J, const auto& gp) {
      const auto geo          = underlying().localView().element().geometry();
      const double intElement = geo.integrationElement(gp.position()) * gp.weight();
      K.template block<FE::myDim, FE::myDim>(I * FE::myDim, J * FE::myDim) += kgIJ * intElement;
    };
  }
  /**
   * \brief Get a lambda function that evaluates the material part of the stiffness matrix (Ke + Ku) for a given strain,
   * integration point and its index.
   *
   * \tparam ST The scalar type for the material and strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \param K The matrix to store the calculated result.
   * \return A lambda function that evaluates the stiffness matrix for a given strain, integration point and its index.
   */
  template <typename ST>
  auto materialStiffnessMatrixFunction(const Requirement& par, typename Traits::template MatrixType<ST>& K,
                                       const VectorXOptRef<ST>& dx = std::nullopt) const {
    return
        [&](const auto& C, const BopType<ST>& bopI, const BopType<ST>& bopJ, const int I, const int J, const auto& gp) {
          const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
          K.template block<myDim, myDim>(I * myDim, J * myDim) += (bopI.transpose() * C * bopJ) * intElement;
        };
  }

  /**
   * \brief Get a lambda function that evaluates the internal force vector for a given strain, integration point and its
   * index.
   *
   * \tparam ST The scalar type for the material and strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \param force The vector to store the calculated result.
   * \return A lambda function that evaluates the intenral force vector for a given strain, integration point and its
   * index.
   */
  template <typename ST>
  auto internalForcesFunction(const Requirement& par, typename Traits::template VectorType<ST>& force,
                              const VectorXOptRef<ST>& dx = std::nullopt) const {
    return [&](const StrainType<ST>& stresses, const BopType<ST>& bopI, const int I, const auto& gp) {
      const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
      force.template segment<myDim>(myDim * I) += bopI.transpose() * stresses * intElement;
    };
  }

  /**
   * \brief Get a lambda function that evaluates the internal energy at a given integration point and its index.
   *
   * \tparam ST The scalar type for the material and strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return A lambda function that returns the intenral energy at a given integration point and its index.
   */
  template <typename ST>
  auto energyFunction(const Requirement& par, const VectorXOptRef<ST>& dx = std::nullopt) const {
    return [&]() -> ST {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      ST energy            = 0.0;
      const auto eps       = strainFunction(par, dx);
      const auto pFunction = pressureFunction(par, dx);
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const auto EVoigt = eps.evaluate(gpIndex, on(gridElement));
        const auto p      = pFunction.evaluate(gpIndex, on(gridElement)).eval();

        auto e = constitutiveDriver<ST>().storedEnergy(EVoigt, p[0]);
        energy += e * geo_->integrationElement(gp.position()) * gp.weight();
      }
      return energy;
    };
  }

protected:
  /**
   * \brief Calculate the matrix associated with the given Requirement.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param par The Requirement object specifying the requirements for the calculation.
   * \param K The matrix to store the calculated result.
   */
  template <typename ScalarType>
  void calculateMatrixImpl(const Requirement& par, const MatrixAffordance& affordance,
                           typename Traits::template MatrixType<> K,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    if constexpr (Concepts::AutodiffScalar<ScalarType>) {
      static_assert(Dune::AlwaysFalse<ScalarType>::value,
                    "DisplacementPressure element does not support matrix calculations for autodiff scalars");
    }
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto pFunction = pressureFunction(par, dx);

    auto pFE         = underlying().localView().tree().child(Dune::Indices::_1);
    auto plocalBasis = pFE.finiteElement().localBasis();
    std::vector<Dune::FieldVector<double, 1>> Np;

    const auto eps        = strainFunction(par, dx);
    const auto kMFunction = materialStiffnessMatrixFunction<ScalarType>(par, K, dx);
    const auto kGFunction = geometricStiffnessMatrixFunction<ScalarType>(par, K, dx);
    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const auto EVoigt = (eps.evaluate(gpIndex, on(gridElement))).eval();
      const auto p      = pFunction.evaluate(gpIndex, on(gridElement)).eval();
      plocalBasis.evaluateFunction(gp.position(), Np);
      const auto intElement = geo_->integrationElement(gp.position()) * gp.weight();

      auto [dUdE, S] = constitutiveDriver<ScalarType>().firstDerivatives(EVoigt, p[0]);
      auto C         = constitutiveDriver<ScalarType>().secondDerivative(EVoigt, p[0]);

      // Loop over u-u
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        for (size_t j = 0; j < numberOfNodes_; ++j) {
          const auto bopJ = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), on(gridElement));
          const auto kgIJ = eps.evaluateDerivative(gpIndex, wrt(coeff(i, j)), along(S), on(gridElement));
          kMFunction(C, bopI, bopJ, i, j, gp);
          kGFunction(kgIJ, i, j, gp);
        }
      }

      // Loop over u-p
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        for (auto j : Dune::range(pFE.size())) {
          auto jIdx         = pFE.localIndex(j);
          const auto upTerm = (bopI.transpose() * dUdE * Np[j][0]).eval();
          K.template block<myDim, 1>(i * myDim, jIdx) += upTerm * intElement;
        }
      }

      // Loop over p-u
      for (auto i : Dune::range(pFE.size())) {
        auto iIdx = pFE.localIndex(i);
        for (size_t j = 0; j < numberOfNodes_; ++j) {
          const auto bopJ   = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), on(gridElement));
          const auto upTerm = (bopJ.transpose() * dUdE * Np[i][0]).transpose().eval();
          K.template block<1, myDim>(iIdx, j * myDim) += upTerm * intElement;
        }
      }

      // Loop over p-p
      for (auto i : Dune::range(pFE.size())) {
        auto iIdx = pFE.localIndex(i);
        for (auto j : Dune::range(pFE.size())) {
          auto jIdx   = pFE.localIndex(j);
          auto ppTerm = Np[i][0] * (1.0 / constitutiveDriver_.kappa()) * Np[j][0];
          K(iIdx, jIdx) -= ppTerm * intElement;
        }
      }
    }
  }

  template <typename ScalarType>
  auto calculateScalarImpl(const Requirement& par, ScalarAffordance affordance,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const -> ScalarType {
    if constexpr (not Concepts::AutodiffScalar<ScalarType>) {
      return energyFunction(par, dx)();
    } else {
      static_assert(Dune::AlwaysFalse<ScalarType>::value,
                    "DisplacementPressure element does not support scalar calculations for autodiff scalars");
    }
  }

  template <typename ScalarType>
  void calculateVectorImpl(const Requirement& par, VectorAffordance affordance,
                           typename Traits::template VectorType<ScalarType> force,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    if constexpr (Concepts::AutodiffScalar<ScalarType>) {
      static_assert(Dune::AlwaysFalse<ScalarType>::value,
                    "DisplacementPressure element does not support vector calculations for autodiff scalars");
    }
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto eps       = strainFunction(par, dx);
    const auto pFunction = pressureFunction(par, dx);

    auto pFE         = underlying().localView().tree().child(Dune::Indices::_1);
    auto plocalBasis = pFE.finiteElement().localBasis();
    std::vector<Dune::FieldVector<double, 1>> Np;

    const auto fIntFunction = internalForcesFunction<ScalarType>(par, force, dx);

    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const auto EVoigt = (eps.evaluate(gpIndex, on(gridElement))).eval();
      const auto p      = pFunction.evaluate(gpIndex, on(gridElement)).eval();
      plocalBasis.evaluateFunction(gp.position(), Np);

      auto [dUdE, S] = constitutiveDriver<ScalarType>().firstDerivatives(EVoigt, p[0]);
      auto Uhat      = constitutiveDriver<ScalarType>().Uhat(EVoigt);

      // Loop over u
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        fIntFunction(S, bopI, i, gp);
      }

      // Loop over p
      for (auto i : Dune::range(pFE.size())) {
        auto iIdx  = pFE.localIndex(i);
        auto pTerm = Np[i][0] * (Uhat - p[0] / constitutiveDriver_.kappa());
        force(iIdx) += pTerm * geo_->integrationElement(gp.position()) * gp.weight();
      }
    }
  }
};

/**
 * \brief A helper function to create a displacement-pressure pre finite element.
 * \tparam MAT Type of the material.
 * \param mat The material model.
 * \return A displacement-pressure pre finite element.
 */
template <Concepts::Material MAT>
auto displacementPressure(const MAT& mat) {
  DisplacementPressurePre<MAT> pre(mat);
  return pre;
}

} // namespace Ikarus

#else
  #error DisplacementPressure depends on dune-localfefunctions, which is not included
#endif

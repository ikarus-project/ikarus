// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <concepts>
#include <iosfwd>
#include <optional>

#include <dune/common/classname.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/expressions/greenLagrangeStrains.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feBases/powerBasisFE.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/finiteElements/mechanics/materials.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename Basis, typename Material>
  class NonLinearElastic : public PowerBasisFE<typename Basis::FlatBasis> {
  public:
    using FlatBasis                  = typename Basis::FlatBasis;
    using BasePowerFE                = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType          = FErequirements<Eigen::VectorXd>;
    using ResultRequirementsType     = ResultRequirements<Eigen::VectorXd>;
    using LocalView                  = typename FlatBasis::LocalView;
    using Geometry                   = typename LocalView::Element::Geometry;
    using GridView                   = typename FlatBasis::GridView;
    using Traits                     = TraitsFromLocalView<LocalView>;
    static constexpr int myDim       = Traits::mydim;
    static constexpr auto strainType = StrainTags::greenLagrangian;

    template <typename VolumeLoad = std::nullptr_t, typename NeumannBoundaryLoad = std::nullptr_t>
    requires(Std::is_pointer<VolumeLoad>and Std::is_pointer<NeumannBoundaryLoad>)
        NonLinearElastic(Basis& globalBasis, const typename LocalView::Element& element, const Material& p_mat,
                              const BoundaryPatch<GridView>* neumannBoundary   = nullptr,
                              const NeumannBoundaryLoad& p_neumannBoundaryLoad = nullptr,
                              const VolumeLoad& p_volumeLoad                   = nullptr)
        : BasePowerFE(globalBasis.flat(), element),
          volumeLoad{Std::returnReferenceOrNulloptIfObjectIsNullPtr(p_volumeLoad)},
          neumannBoundaryLoad_{Std::returnReferenceOrNulloptIfObjectIsNullPtr(p_neumannBoundaryLoad)},
          neumannBoundary_{neumannBoundary},
          mat{p_mat} {
      this->localView().bind(element);
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes     = fe.size();
      dispAtNodes.resize(fe.size());
      const int order = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      localBasis      = Dune::CachedLocalBasis(this->localView().tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(this->localView().element().type(), order),
                      Dune::bindDerivatives(0, 1));
      assert(((not neumannBoundary_ and not neumannBoundaryLoad_) or (neumannBoundary_ and neumannBoundaryLoad_))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

  public:
    template <typename ScalarType>
    auto getDisplacementFunction(const FERequirementType& par, const Eigen::VectorX<ScalarType>& dx) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);

      Dune::BlockVector<Dune::RealTuple<ScalarType, Traits::dimension>> disp(dispAtNodes.size());
      for (auto i = 0U; i < disp.size(); ++i)
        for (auto k2 = 0U; k2 < myDim; ++k2)
          disp[i][k2]
              = dx[i * myDim + k2] + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];

      auto geo = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          this->localView().element().geometry());
      Dune::StandardLocalFunction uFunction(localBasis, disp, geo);
      return uFunction;
    }

    template <typename ScalarType = double>
    auto getStrainFunction(const FERequirementType& par, const Eigen::VectorX<ScalarType>& dx) const {
      return greenLagrangeStrains(getDisplacementFunction(par, dx));
    }

    template <bool voigt = true>
    auto getMaterialTangent(const auto& strain) const {
      return mat.template tangentModuli<strainType, voigt>(strain);
    }

    template <typename ScalarType, int strainDim>
    auto getInternalEnergy(const Eigen::Vector<ScalarType, strainDim>& strain) const {
      if constexpr (std::is_same_v<ScalarType, double>)
        return mat.template storedEnergy<strainType>(strain);
      else {
        decltype(auto) matAD = mat.template rebind<ScalarType>();
        return matAD.template storedEnergy<strainType>(strain);
      }
    }

    template <bool voigt = true>
    auto getStress(const auto& strain) const {
      return mat.template stresses<strainType, voigt>(strain);
    }

    template <typename ScalarType = double>
    ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(this->localView().size());
      dx.setZero();
      return calculateScalarImpl(par, dx);
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& K) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      Eigen::VectorXd dx(this->localView().size());
      dx.setZero();
      const auto eps = getStrainFunction(par, dx);
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

    void calculateAt(const ResultRequirementsType& req, const Eigen::Vector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::Indices;
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      Eigen::VectorXd dx(this->localView().size());
      dx.setZero();

      const auto eps    = getStrainFunction(req.getFERequirements(), dx);
      const auto gp     = toDune(local);
      const auto EVoigt = (eps.evaluate(gp, on(gridElement))).eval();

      typename ResultTypeMap<double>::ResultArray resultVector;
      if (req.isResultRequested(ResultType::secondPiolaKirchhoffStress)) {
        const auto stress = getStress<false>(EVoigt);
        resultVector.resizeLike(stress);
        resultVector = stress;
        result.insertOrAssignResult(ResultType::secondPiolaKirchhoffStress, resultVector);
      }
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& force) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      Eigen::VectorXd dx(this->localView().size());
      dx.setZero();
      const auto eps = getStrainFunction(par, dx);
      const auto geo = this->localView().element().geometry();

      // Internal forces
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        const auto EVoigt       = (eps.evaluate(gpIndex, on(gridElement))).eval();
        const auto C            = getMaterialTangent(EVoigt);
        const auto stresses     = getStress(EVoigt);
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          force.template segment<myDim>(myDim * i) += bopI.transpose() * stresses * intElement;
        }
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        const auto u = getDisplacementFunction(par, dx);
        for (const auto& [gpIndex, gp] : u.viewOverIntegrationPoints()) {
          const double intElement                            = geo.integrationElement(gp.position()) * gp.weight();
          const Eigen::Vector<double, Traits::worlddim> fExt = (*volumeLoad)(toEigen(gp.position()), lambda);
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<myDim>(myDim * i) -= udCi * fExt * intElement;
          }
        }
      }

      // External forces, boundary forces, i.e. at the Neumann boundary
      if (not neumannBoundary_ and not neumannBoundaryLoad_) return;

      const auto u        = getDisplacementFunction(par, dx);
      const auto& element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
        if (not neumannBoundary_->contains(intersection)) continue;

        // Integration rule along the boundary
        const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), u.order());

        for (const auto& curQuad : quadLine) {
          const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double intElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function wrt the i-th coefficient
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(quadPos, wrt(coeff(i)));

            // Value of the Neumann data at the current position
            const auto neumannValue
                = (*neumannBoundaryLoad_)(toEigen(intersection.geometry().global(curQuad.position())), lambda);
            force.template segment<myDim>(myDim * i) -= udCi * neumannValue * curQuad.weight() * intElement;
          }
        }
      }
    }

    Dune::CachedLocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::optional<std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                                        const double&)>>
        volumeLoad;
    std::optional<std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                                        const double&)>>
        neumannBoundaryLoad_;
    const BoundaryPatch<GridView>* neumannBoundary_;
    mutable Dune::BlockVector<Dune::RealTuple<double, Traits::dimension>> dispAtNodes;
    Material mat;
    size_t numberOfNodes{0};

  protected:
    template <typename ScalarType = double>
    ScalarType calculateScalarImpl(const FERequirementType& par, const Eigen::VectorX<ScalarType>& dx) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = getDisplacementFunction(par, dx);
      const auto eps       = getStrainFunction(par, dx);
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
          const Eigen::Vector<double, Traits::worlddim> fExt = (*volumeLoad)(toEigen(gp.position()), lambda);
          energy -= u.dot(fExt) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary_ and not neumannBoundaryLoad_) return energy;

      const auto& element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
        if (not neumannBoundary_ or not neumannBoundary_->contains(intersection)) continue;

        const auto& quadLine
            = Dune::QuadratureRules<double, Traits::mydim - 1>::rule(intersection.type(), uFunction.order());

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, Traits::mydim>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double intElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto u = uFunction.evaluate(quadPos);

          // Value of the Neumann data at the current position
          const auto neumannValue
              = (*neumannBoundaryLoad_)(toEigen(intersection.geometry().global(curQuad.position())), lambda);
          energy -= neumannValue.dot(u) * curQuad.weight() * intElement;
        }
      }
      return energy;
    }

  public:
    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      const auto& d      = req.getGlobalSolution(Ikarus::FESolutions::displacement);
      const auto& lambda = req.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Dune::BlockVector<RealTuple<double, Traits::dimension>> disp(fe.size());

      for (auto i = 0U; i < fe.size(); ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      const auto geo = this->localView().element().geometry();

      Dune::StandardLocalFunction uFunction(localBasis, disp, std::make_shared<const Geometry>(geo));
      const auto H      = uFunction.evaluateDerivative(local, Dune::wrt(spatialAll), Dune::on(gridElement));
      const auto E      = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
      const auto EVoigt = toVoigt(E);
      auto PK2          = mat.template stresses<StrainTags::greenLagrangian>(EVoigt);

      typename ResultTypeMap<double>::ResultArray resultVector;
      if (req.isResultRequested(ResultType::PK2Stress)) {
        resultVector.resizeLike(PK2);
        resultVector = PK2;
        result.insertOrAssignResult(ResultType::PK2Stress, resultVector);
      }
    }

  private:
    Dune::CachedLocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        volumeLoad;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        neumannBoundaryLoad;
    const BoundaryPatch<GridView>* neumannBoundary;
    Material mat;
  };
}  // namespace Ikarus

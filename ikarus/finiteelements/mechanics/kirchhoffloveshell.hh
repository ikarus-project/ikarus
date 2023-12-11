// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

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

#include <ikarus/finiteelements/febases/powerbasisfe.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/fetraits.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/eigendunetransformations.hh>
#include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

  template <class ScalarType>
  ScalarType energyHelper(const Eigen::Vector<ScalarType, 3>& epsV, const auto& Aconv, double E, double nu) {
    const double lambda   = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu       = E / (2.0 * (1.0 + nu));
    const double lambdbar = 2.0 * lambda * mu / (lambda + 2.0 * mu);
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>> moduli;
    const auto AconvT = TensorCast(Aconv, std::array<Eigen::Index, 2>({3, 3}));
    moduli            = lambdbar * dyadic(AconvT, AconvT).eval() + 2 * mu * symmetricFourthOrder<double>(Aconv, Aconv);

    auto C   = toVoigt(moduli);
    auto C33 = C({0, 1, 5}, {0, 1, 5}).eval();

    return 0.5 * epsV.dot(C33 * epsV);
  }

  template <typename Basis_, typename FERequirements_ = FErequirements<>, bool useEigenRef = false>
  class KirchhoffLoveShell : public PowerBasisFE<typename Basis_::FlatBasis> {
  public:
    using Basis                   = Basis_;
    using FlatBasis               = typename Basis::FlatBasis;
    using BasePowerFE             = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType       = FERequirements_;
    using ResultRequirementsType  = ResultRequirements<FERequirementType>;
    using LocalView               = typename FlatBasis::LocalView;
    using Element                 = typename LocalView::Element;
    using Geometry                = typename Element::Geometry;
    using GridView                = typename FlatBasis::GridView;
    using Traits                  = TraitsFromLocalView<LocalView, useEigenRef>;
    static constexpr int myDim    = Traits::mydim;
    static constexpr int worlddim = Traits::worlddim;

    template <typename VolumeLoad = LoadDefault, typename NeumannBoundaryLoad = LoadDefault>
    KirchhoffLoveShell(const Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu,
                       double thickness, VolumeLoad p_volumeLoad = {},
                       const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
                       NeumannBoundaryLoad p_neumannBoundaryLoad        = {})
        : BasePowerFE(globalBasis.flat(), element),
          neumannBoundary{p_neumannBoundary},
          emod_{emod},
          nu_{nu},
          thickness_{thickness} {
      this->localView().bind(element);
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes     = fe.size();
      dispAtNodes.resize(fe.size());
      order      = 2 * (fe.localBasis().order());
      localBasis = Dune::CachedLocalBasis(fe.localBasis());
      if constexpr (requires { this->localView().element().impl().getQuadratureRule(order); })
        if (this->localView().element().impl().isTrimmed())
          localBasis.bind(this->localView().element().impl().getQuadratureRule(order), Dune::bindDerivatives(0, 1, 2));
        else
          localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                          Dune::bindDerivatives(0, 1, 2));
      else
        localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                        Dune::bindDerivatives(0, 1, 2));

      if constexpr (!std::is_same_v<VolumeLoad, LoadDefault>) volumeLoad = p_volumeLoad;
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, LoadDefault>) neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

  public:
    template <typename ScalarType = double>
    auto getDisplacementFunction(const FERequirementType& par,
                                 const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);

      Dune::BlockVector<Dune::RealTuple<ScalarType, Traits::worlddim>> disp(dispAtNodes.size());
      if (dx)
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            disp[i][k2] = dx.value()[i * worlddim + k2]
                          + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      else
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            disp[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];

      auto geo = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          this->localView().element().geometry());
      Dune::StandardLocalFunction uFunction(localBasis, disp, geo);
      return std::make_pair(uFunction, disp);
    }

    inline double calculateScalar(const FERequirementType& par) const { return calculateScalarImpl<double>(par); }

    void calculateAt([[maybe_unused]] const ResultRequirementsType& req,
                     [[maybe_unused]] const Dune::FieldVector<double, Traits::mydim>& local,
                     [[maybe_unused]] ResultTypeMap<double>& result) const {
      DUNE_THROW(Dune::NotImplemented, "No results are implemented");
    }

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
    mutable Dune::BlockVector<Dune::RealTuple<double, Traits::dimension>> dispAtNodes;
    double emod_;
    double nu_;
    double thickness_;
    size_t numberOfNodes{0};
    int order{};

  protected:
    template <typename ScalarType>
    auto calculateScalarImpl(const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                                           = std::nullopt) const -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [uFunction, uNodes] = getDisplacementFunction(par, dx);
      const auto& lambda             = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto geo                 = this->localView().element().geometry();
      ScalarType energy              = 0.0;
      const auto uasMatrix           = Dune::viewAsEigenMatrixAsDynFixed(uNodes);

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto [X, Jd, Hd]                      = geo.impl().zeroFirstAndSecondDerivativeOfPosition(gp.position());
        const auto J                                = toEigen(Jd);
        const auto H                                = toEigen(Hd);
        const Eigen::Matrix<double, 2, 2> A         = J * J.transpose();
        const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
            uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, Dune::on(DerivativeDirections::referenceElement))));
        const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();

        const auto& Ndd                     = localBasis.evaluateSecondDerivatives(gpIndex);
        const auto h                        = H + Ndd.transpose().template cast<ScalarType>() * uasMatrix;
        const Eigen::Vector3<ScalarType> a3 = (j.row(0).cross(j.row(1))).normalized();
        Eigen::Vector<ScalarType, 3> bV     = h * a3;
        bV(2) *= 2;  // Voigt notation requires the two here

        Eigen::Matrix<double, 3, 3> G;
        G.setZero();
        G.block<2, 2>(0, 0)                    = A;
        G(2, 2)                                = 1;
        const Eigen::Matrix<double, 3, 3> GInv = G.inverse();

        const auto epsV                 = toVoigt((0.5 * (j * j.transpose() - A)).eval()).eval();
        const auto BV                   = toVoigt(toEigen(geo.impl().secondFundamentalForm(gp.position())));
        const auto kappaV               = (BV - bV).eval();
        const ScalarType membraneEnergy = thickness_ * energyHelper(epsV, GInv, emod_, nu_);
        const ScalarType bendingEnergy  = Dune::power(thickness_, 3) / 12.0 * energyHelper(kappaV, GInv, emod_, nu_);
        energy += (membraneEnergy + bendingEnergy) * geo.integrationElement(gp.position()) * gp.weight();
      }

      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
          const auto u                                       = uFunction.evaluate(gpIndex);
          const Eigen::Vector<double, Traits::worlddim> fExt = volumeLoad(toEigen(geo.global(gp.position())), lambda);
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
          const auto neumannValue
              = neumannBoundaryLoad(toEigen(intersection.geometry().global(curQuad.position())), lambda);
          energy -= neumannValue.dot(u) * curQuad.weight() * intElement;
        }
      }

      return energy;
    }
  };
}  // namespace Ikarus

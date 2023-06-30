// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>
#include <iosfwd>
#include <numbers>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/impl/projectionBasedLocalFunction.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/localfefunctions/manifolds/realTuple.hh>
#include <dune/localfefunctions/manifolds/unitVector.hh>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/utils/defaultFunctions.hh>
#include <dune/fufem/boundarypatch.hh>
#include "fesettings.hh"

namespace Ikarus {

  template <typename Basis_, bool useEigenRef_ = false>
  class NonLinearRMshell {
   public:
    using Basis = Basis_;

    template<typename ScalarType>
    struct KinematicVariables {
      // current configuration
      Eigen::Matrix<ScalarType, 3, 2> a1anda2;  // first and second tangent base vector on the current geometry
      Eigen::Matrix<ScalarType, 3, 2>
          td1Andtd2;  // partial derivative of the director on the current geometry in first and second direction
      Eigen::Matrix<ScalarType, 3, 2>
          ud1andud2;  // partial derivative of the displacement function in first and second direction
      Eigen::Vector<ScalarType, 3> t;  // unit director on the current geometry

      // reference configuration
      Eigen::Matrix<double, 3, 2> A1andA2;  // first and second tangent base vector on the reference geometry
      Eigen::Matrix<double, 3, 2>
          t0d1Andt0d2;  // partial derivative of the director on the reference geometry in first and second direction
      Eigen::Vector<double, 3> t0;  // unit director on the reference geometry

      auto A1() const { return A1andA2.col(0); }
      auto A2() const { return A1andA2.col(1); }
      auto a1() const { return a1anda2.col(0); }
      auto a2() const { return a1anda2.col(1); }
      auto ud1() const { return ud1andud2.col(0); }
      auto ud2() const { return ud1andud2.col(1); }

      auto td1() const { return td1Andtd2.col(0); }
      auto td2() const { return td1Andtd2.col(1); }

      auto t0d1() const { return t0d1Andt0d2.col(0); }
      auto t0d2() const { return t0d1Andt0d2.col(1); }
    };

  public:
    static_assert(Basis::UntouchedBasis::PreBasis::template SubPreBasis<0>::Node::degree() == 3,
                  "The midsurface positions need a power 3D basis");
    static_assert(Basis::UntouchedBasis::PreBasis::template SubPreBasis<1>::Node::degree() == 2,
                  "The director corrections need a power 2D basis");

    static constexpr int midSurfaceDim         = 3;
    static constexpr int directorDim           = 3;
    static constexpr int directorCorrectionDim = directorDim - 1;
    using MidSufaceVector                      = Dune::BlockVector<Dune::RealTuple<double, midSurfaceDim>>;
    using DirectorVector                       = Dune::BlockVector<Dune::UnitVector<double, directorDim>>;
    using MultiTypeVector                      = Dune::TupleVector<MidSufaceVector, DirectorVector>;

    using FERequirementType      = FErequirements<MultiTypeVector>;
    using ResultRequirementsType = ResultRequirements<MultiTypeVector>;
    using LocalViewBlocked       = typename Basis::UntouchedBasis::LocalView;
    using LocalViewFlat          = typename Basis::FlatBasis::LocalView;
    using Geometry               = typename LocalViewFlat::Element::Geometry;
    using GridView          = typename Basis::FlatBasis::GridView;
    using Traits = TraitsFromLocalView<LocalViewFlat, false>;
    static constexpr int useEigenRef = useEigenRef_;

    static constexpr int myDim = Traits::mydim;
    static constexpr int worlddim = Traits::worlddim;

    template <typename VolumeLoad = LoadDefault, typename NeumannBoundaryLoad = LoadDefault>
    NonLinearRMshell(const Basis &globalBasis, const typename LocalViewBlocked::Element &element,
                     const FESettings &p_feSettings,  const MultiTypeVector &p_x0,
        VolumeLoad p_volumeLoad = {}, const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
        NeumannBoundaryLoad p_neumannBoundaryLoad = {},
                     const double kappa_ = 1.0)
        : localViewFlat_{globalBasis.flat().localView()},
          localViewBlocked_{globalBasis.untouched().localView()},
         fESettings{std::move(p_feSettings)},
          x0{p_x0},
          neumannBoundary(p_neumannBoundary)
         {
           if constexpr (!std::is_same_v<VolumeLoad, LoadDefault>) volumeLoad = p_volumeLoad;
           if constexpr (!std::is_same_v<NeumannBoundaryLoad, LoadDefault>) neumannBoundaryLoad = p_neumannBoundaryLoad;
      localViewFlat_.bind(element);
      localViewBlocked_.bind(element);
      geo_ = std::make_shared<const Geometry>(localViewFlat_.element().geometry());

      const auto &child0 = localViewBlocked_.tree().child(Dune::Indices::_0, 0);
      const auto &fe0    = child0.finiteElement();
      const auto &child1 = localViewBlocked_.tree().child(Dune::Indices::_1, 0);
      const auto &fe1    = child1.finiteElement();

      numNodeMidSurface = fe0.size();
      numNodeDirector   = fe1.size();

      order = std::max(2 * localViewFlat_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order(),
                       2 * localViewFlat_.tree().child(Dune::Indices::_1, 0).finiteElement().localBasis().order());
      localBasisMidSurface
          = Dune::CachedLocalBasis(localViewFlat_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis());
      localBasisMidSurface.bind(
          Dune::QuadratureRules<double, Traits::mydim>::rule(localViewFlat_.element().type(), order),
          Dune::bindDerivatives(0, 1));
      localBasisDirector
          = Dune::CachedLocalBasis(localViewFlat_.tree().child(Dune::Indices::_1, 0).finiteElement().localBasis());
      localBasisDirector.bind(
          Dune::QuadratureRules<double, Traits::mydim>::rule(localViewFlat_.element().type(), order),
          Dune::bindDerivatives(0, 1));
      const auto &Emodul = fESettings.request<double>("youngs_modulus");
      const auto &nu = fESettings.request<double>("poissons_ratio");
      thickness_ = fESettings.request<double>("thickness");

      CMat_.setZero();
      // membrane
      const double fac1 = thickness_ * Emodul / (1 - nu * nu);
      CMat_(0, 0) = CMat_(1, 1) = fac1;
      CMat_(2, 2)               = fac1 * (1 - nu) * 0.5;
      CMat_(1, 0) = CMat_(0, 1) = fac1 * nu;

      // bending
      const double fac2 = thickness_ * thickness_ * thickness_ / 12 * Emodul / (1 - nu * nu);
      CMat_(3, 3) = CMat_(4, 4) = fac2;
      CMat_(5, 5)               = fac2 * (1 - nu) * 0.5;
      CMat_(3, 4) = CMat_(4, 3) = fac2 * nu;

      // transverse shear
      const double fac3 = kappa_ * thickness_ * Emodul * 0.5 / (1 + nu);
      CMat_(6, 6) = CMat_(7, 7) = fac3;

      if constexpr (!std::is_same_v<VolumeLoad, LoadDefault>) volumeLoad = p_volumeLoad;
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, LoadDefault>) neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
                 && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }
    FESettings fESettings;

    int ndofEmbedded()const
    {
      using namespace Dune::Indices;
      return localViewFlat_.tree().child(_0, 0).finiteElement().size()*3+ localViewFlat_.tree().child(_1, 0).finiteElement().size()*3;
    }

    using GlobalIndex = typename LocalViewFlat::MultiIndex;
    void globalFlatIndices(std::vector<GlobalIndex> &globalIndices) const {
      using namespace Dune::Indices;
      const auto &fe = localViewFlat_.tree().child(_0, 0).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < midSurfaceDim; ++j) {
          globalIndices.push_back(localViewFlat_.index((localViewFlat_.tree().child(_0, j).localIndex(i))));
        }
      }
      const auto &fe2 = localViewFlat_.tree().child(_1, 0).finiteElement();
      for (size_t i = 0; i < fe2.size(); ++i) {
        for (int j = 0; j < directorCorrectionDim; ++j) {
          globalIndices.push_back(localViewFlat_.index((localViewFlat_.tree().child(_1, j).localIndex(i))));
        }
      }
    }

    template<typename ScalarType>
    static Eigen::Vector<ScalarType, 8> calculateGreenLagrangianStrains(const KinematicVariables<ScalarType> &kin) {
      Eigen::Vector<ScalarType, 8> egl;
      // membrane
      egl << kin.A1().dot(kin.ud1()) + 0.5 * kin.ud1().squaredNorm(),
          kin.A2().dot(kin.ud2()) + 0.5 * kin.ud2().squaredNorm(), kin.a1().dot(kin.a2()),
          // bending
          kin.a1().dot(kin.td1()) - kin.A1().dot(kin.t0d1()), kin.a2().dot(kin.td2()) - kin.A2().dot(kin.t0d2()),
          kin.a1().dot(kin.td2()) + kin.a2().dot(kin.td1()) - kin.A1().dot(kin.t0d2()) - kin.A2().dot(kin.t0d1()),
          // trans shear
          kin.a1().dot(kin.t) - kin.A1().dot(kin.t0), kin.a2().dot(kin.t) - kin.A2().dot(kin.t0);

      return egl;
    }

    template <typename ScalarType>
    static Eigen::Vector<ScalarType, 8> calculateStressResultants(const Eigen::Matrix<double, 8, 8> &Cred,
                                                              const Eigen::Vector<ScalarType, 8> &Egl) {
      Eigen::Vector<ScalarType, 8> S;
      S << Cred(0, 0) * Egl(0) + Cred(0, 1) * Egl(1), Cred(1, 0) * Egl(0) + Cred(1, 1) * Egl(1), Cred(2, 2) * Egl(2),
          Cred(3, 3) * Egl(3) + Cred(3, 4) * Egl(4), Cred(4, 3) * Egl(3) + Cred(4, 4) * Egl(4), Cred(5, 5) * Egl(5),
          Cred(6, 6) * Egl(6), Cred(7, 7) * Egl(7);

      return S;
    }

    auto size() const { return localViewFlat_.size(); }
    auto& localViewBlocked() const { return localViewBlocked_; }

    auto getReferenceLocalConfigurations() const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const int nDofs0 = localViewBlocked_.tree().child(_0, 0).finiteElement().size();
      const int nDofs1 = localViewBlocked_.tree().child(_1, 0).finiteElement().size();

//      std::vector<Dune::RealTuple<double, 3>> localConfiguration0(nDofs0);
      std::vector<Dune::UnitVector<double, 3>> localConfiguration1(nDofs1);

      for (int i = 0; i < nDofs0 + nDofs1; i++) {
        int localIndexI = 0;
        if (i < nDofs0) {
          auto &node  = localViewBlocked_.tree().child(_0, 0);
          localIndexI = node.localIndex(i);
        } else {
          auto &node  = localViewBlocked_.tree().child(_1, 0);
          localIndexI = node.localIndex(i - nDofs0);
        }
        auto multiIndex = localViewBlocked_.index(localIndexI);

        // The CompositeBasis number is contained in multiIndex[0]
        // multiIndex[1] contains the actual index
//        if (multiIndex[0] == 0)
//          localConfiguration0[i] = x0[_0][multiIndex[1]];
//        else
          if (multiIndex[0] == 1)
          localConfiguration1[i - nDofs0] = x0[_1][multiIndex[1]];
      }
      return localConfiguration1;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 8, 3> boperatorMidSurface(const KinematicVariables<ScalarType> &kin, int integrationPointIndex,
                                                        int coeffIndex, const auto &midSurfaceFunction) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda2
          = midSurfaceFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll, coeff(_0, coeffIndex)));
      const auto &diffa1 = diffa1Anda2[0];
      const auto &diffa2 = diffa1Anda2[1];
      Eigen::Matrix<ScalarType, 8, 3> bop;
      bop.setZero();
      bop.row(0) = kin.a1().transpose() * diffa1;  // membrane_{,disp}
      bop.row(1) = kin.a2().transpose() * diffa2;
      bop.row(2) = kin.a2().transpose() * diffa1 + kin.a1().transpose() * diffa2;

      bop.row(3) = kin.td1().transpose() * diffa1;  // bending_{,disp}
      bop.row(4) = kin.td2().transpose() * diffa2;
      bop.row(5) = kin.td2().transpose() * diffa1 + kin.td1().transpose() * diffa2;

      bop.row(6) = kin.t.transpose() * diffa1;  // trans_shear_{,disp}
      bop.row(7) = kin.t.transpose() * diffa1;
      return bop;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 8, 2> boperatorDirector(const KinematicVariables<ScalarType> &kin, int integrationPointIndex,
                                                      int coeffIndex, const auto &directorFunction) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const std::array<Eigen::Matrix<ScalarType, 3, 2>, 2> diffdt
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll, coeff(_1, coeffIndex)));
      const Eigen::Matrix<ScalarType, 3, 2> difft
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(coeff(_1, coeffIndex)));

      const auto &diffdtd1 = diffdt[0];
      const auto &diffdtd2 = diffdt[1];

      Eigen::Matrix<ScalarType, 8, 2> bop;
      bop.setZero();
      bop.row(3) = kin.a1().transpose() * diffdtd1;  // bending_{,disp}
      bop.row(4) = kin.a2().transpose() * diffdtd2;
      bop.row(5) = kin.a2().transpose() * diffdtd1 + kin.a1().transpose() * diffdtd2;

      bop.row(6) = kin.a1().transpose() * difft;  // trans_shear_{,disp}
      bop.row(7) = kin.a2().transpose() * difft;

      return bop;
    }

    template<typename ScalarType>
    auto createFunctions(const FERequirementType &par,
                         const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
      using namespace Dune::Indices;

      const auto &mNodal = par.getGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector)[_0];
      const auto &dNodal = par.getGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector)[_1];
      //      const auto &lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      const auto &child0 = localViewBlocked_.tree().child(_0, 0);
      const auto &fe0    = child0.finiteElement();
      const auto &child1 = localViewBlocked_.tree().child(_1, 0);
      const auto &fe1    = child1.finiteElement();

      Dune::BlockVector<typename std::remove_cvref_t<decltype(dNodal[0])>::template rebind<ScalarType>::other> localDirectorConfiguration(fe1.size());
      Dune::BlockVector<typename std::remove_cvref_t<decltype(mNodal[0])>::template rebind<ScalarType>::other>   displacements(fe0.size());
      if (dx) {
        for (auto i = 0U; i < fe0.size(); ++i) {
          const auto globalIndex = localViewBlocked_.index(localViewBlocked_.tree().child(_0).localIndex(i));
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            displacements[i][k2] = mNodal[globalIndex[1]][k2] +dx.value()[i*worlddim + k2];
        }
        const int nDofs0 = this->localViewBlocked().tree().child(_0, 0).finiteElement().size()*3;

        for (auto i = 0U; i < fe1.size(); ++i) {
          const auto globalIndex = localViewBlocked_.index(localViewBlocked_.tree().child(_1).localIndex(i));
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            localDirectorConfiguration[i][k2] = dNodal[globalIndex[1]][k2]+dx.value()[nDofs0+i*worlddim + k2];
        }
      }else
      {        for (auto i = 0U; i < fe0.size(); ++i) {
          const auto globalIndex = localViewBlocked_.index(localViewBlocked_.tree().child(_0).localIndex(i));
          displacements[i] = mNodal[globalIndex[1]];
        }

        for (auto i = 0U; i < fe1.size(); ++i) {
          const auto globalIndex = localViewBlocked_.index(localViewBlocked_.tree().child(_1).localIndex(i));
          localDirectorConfiguration[i] = dNodal[globalIndex[1]];
        }}

      const auto localRefDirectorConfiguration = getReferenceLocalConfigurations();


      using namespace Dune::DerivativeDirections;

      Dune::StandardLocalFunction displacementFunction(localBasisMidSurface, displacements, geo_, _0);
      Dune::ProjectionBasedLocalFunction directorFunction(localBasisDirector, localDirectorConfiguration, geo_, _1);
      Dune::ProjectionBasedLocalFunction directorReferenceFunction(localBasisDirector, localRefDirectorConfiguration,
                                                                   geo_, _1);

      return std::make_tuple( displacementFunction, directorFunction,
                             directorReferenceFunction);
    }

//    template<typename ScalarType>
//    void calculateMatrixImpl(const FERequirementType &par, typename Traits::template MatrixType<ScalarType> hred,
//                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
//      using namespace Dune::Indices;
//      hred.setZero();
//
//      //      const auto &lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
//
//      using namespace Dune::DerivativeDirections;
//      using namespace Dune::Indices;
//      const auto [midSurfaceFunction, midSurfaceReferenceFunction, displacementFunction, directorFunction,
//                  directorReferenceFunction]
//          = createFunctions(par);
//
//      KinematicVariables kin{};
//
//      Dune::BlockVector<Dune::FieldMatrix<field_type, 8, 3>> bopMidSurface(numNodeMidSurface);
//      Dune::BlockVector<Dune::FieldMatrix<field_type, 8, 2>> bopDirector(numNodeDirector);
//
//      const int midSurfaceDofs = numNodeMidSurface * midSurfaceDim;
//      for (const auto &[gpIndex, gp] : midSurfaceFunction.viewOverIntegrationPoints()) {
//        const auto weight = geo_->integrationElement(gp.position()) * gp.weight();
//
//        kin.t           = directorFunction.evaluate(gpIndex,Dune::on(referenceElement));
//        kin.t0          = directorReferenceFunction.evaluate(gpIndex,Dune::on(referenceElement));
//        kin.a1anda2     = midSurfaceFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//        kin.A1andA2     = midSurfaceReferenceFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//
//        auto Egl = calculateGreenLagrangianStrains(kin);
//        //        const auto StimesWeight = calculateStress(CMat_, Egl) * weight;
//
//        for (int i = 0; i < numNodeMidSurface; ++i) {
//          const auto indexI = midSurfaceDim * i;
//          const auto bopI   = boperatorMidSurface(kin, gpIndex, i, midSurfaceFunction);
//          for (int j = 0; j < numNodeMidSurface; ++j) {
//            const auto indexJ = midSurfaceDim * j;
//
//            const auto bopJ = boperatorMidSurface(kin, gpIndex, j, midSurfaceFunction);
//            //            const auto KgIJ = midSurfaceFunction.evaluateDerivative(gpIndex, Dune::wrt(Dune::coeff(_0, i,
//            //            _0, j)));
//
//            hred.template block<midSurfaceDim, midSurfaceDim>(indexI, indexJ)
//                += bopI.transpose() * CMat_ * bopJ * weight;
//          }
//
//          for (int j = 0; j < numNodeDirector; ++j) {
//            const auto indexJ = midSurfaceDofs + directorCorrectionDim * j;
//            const auto bopJ   = boperatorDirector(kin, gpIndex, j, directorFunction);
//
//            hred.template block<midSurfaceDim, directorCorrectionDim>(indexI, indexJ)
//                += bopI.transpose() * CMat_ * bopJ * weight;
//          }
//        }
//
//        for (int i = 0; i < numNodeDirector; ++i) {
//          const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
//          const auto bopI   = boperatorDirector(kin, gpIndex, i, directorFunction);
//          for (int j = 0; j < numNodeMidSurface; ++j) {
//            const auto indexJ = midSurfaceDim * j;
//
//            const auto bopJ = boperatorMidSurface(kin, gpIndex, j, midSurfaceFunction);
//
//            hred.template block<directorCorrectionDim, midSurfaceDim>(indexI, indexJ)
//                += bopI.transpose() * CMat_ * bopJ * weight;
//          }
//
//          for (int j = 0; j < numNodeDirector; ++j) {
//            const auto indexJ = midSurfaceDofs + directorCorrectionDim * j;
//
//            const auto bopJ = boperatorDirector(kin, gpIndex, j, directorFunction);
//
//            hred.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ)
//                += bopI.transpose() * CMat_ * bopJ * weight;
//          }
//        }
//      }
//    }

    //    [[nodiscard]] int size() const { return localView.size(); }

//    template<typename ScalarType>
//    void calculateVectorImpl(const FERequirementType &par, typename Traits::template VectorType<ScalarType> rieGrad,
//                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
//      using namespace Dune::Indices;
//      rieGrad.setZero();
//
//      const auto &lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
//
//      using namespace Dune::DerivativeDirections;
//      using namespace Dune::Indices;
//      const auto [midSurfaceFunction, midSurfaceReferenceFunction, displacementFunction, directorFunction,
//                  directorReferenceFunction]
//          = createFunctions(par);
//
//      KinematicVariables kin{};
//
//      Dune::BlockVector<Dune::FieldMatrix<field_type, 8, 3>> bopMidSurface(numNodeMidSurface);
//      Dune::BlockVector<Dune::FieldMatrix<field_type, 8, 2>> bopDirector(numNodeDirector);
//      const int midSurfaceDofs = numNodeMidSurface * midSurfaceDim;
//      for (const auto &[gpIndex, gp] : midSurfaceFunction.viewOverIntegrationPoints()) {
//        const auto weight = geo_->integrationElement(gp.position()) * gp.weight();
//
//        kin.t           = directorFunction.evaluate(gpIndex,Dune::on(referenceElement));
//        kin.t0          = directorReferenceFunction.evaluate(gpIndex,Dune::on(referenceElement));
//        kin.a1anda2     = midSurfaceFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//        kin.A1andA2     = midSurfaceReferenceFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
//
//        auto Egl                = calculateGreenLagrangianStrains(kin);
//        const auto StimesWeight = (calculateStressResultants(CMat_, Egl) * weight).eval();
//
//        for (int i = 0; i < numNodeMidSurface; ++i) {
//          const auto indexI = midSurfaceDim * i;
//          const auto bopI   = boperatorMidSurface(kin, gpIndex, i, midSurfaceFunction);
//          rieGrad.template segment<midSurfaceDim>(indexI) += bopI.transpose() * StimesWeight;
//        }
//
//        for (int i = 0; i < numNodeDirector; ++i) {
//          const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
//          const auto bopI   = boperatorDirector(kin, gpIndex, i, directorFunction);
//
//          rieGrad.template segment<directorCorrectionDim>(indexI) += bopI.transpose() * StimesWeight;
//        }
//      }
//
//      // External forces volume forces over the domain
//      if (volumeLoad) {
//        for (const auto& [gpIndex, gp] : displacementFunction.viewOverIntegrationPoints()) {
//          Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(geo_->global(gp.position()), lambda)[0];
//          for (size_t i = 0; i < numNodeMidSurface; ++i) {
//            const auto indexI = midSurfaceDim * i;
//            const auto udCi = displacementFunction.evaluateDerivative(gpIndex, Dune::wrt(coeff(i)));
//            rieGrad.template segment<midSurfaceDim>(indexI)
//                -= udCi * fext * geo_->integrationElement(gp.position()) * gp.weight();
//          }
//        }
//      }
//
//      // External forces, boundary forces, i.e. at the Neumann boundary
//      if (not neumannBoundary) return;
//
//      auto element = localViewFlat_.element();
//      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
//        if (not neumannBoundary->contains(intersection)) continue;
//
//        // Integration rule along the boundary
//        const auto& quadLine = Dune::QuadratureRules<double, 1>::rule(intersection.type(), displacementFunction.order());
//
//        for (const auto& curQuad : quadLine) {
//          const Dune::FieldVector<double, 2>& quadPos = intersection.geometryInInside().global(curQuad.position());
//
//          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());
//
//          // The value of the local function wrt the i-th coef
//          for (size_t i = 0; i < numNodeMidSurface; ++i) {
//            const auto indexI = midSurfaceDim * i;
//            const auto udCi = displacementFunction.evaluateDerivative(quadPos, Dune::wrt(coeff(i)));
//
//            // Value of the Neumann data at the current position
//            auto neumannValue
//                = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda)[0];
//            rieGrad.template segment<midSurfaceDim>(indexI) -= udCi * neumannValue * curQuad.weight() * integrationElement;
//          }
//        }
//      }
//    }

    inline double calculateScalar(const FERequirementType &par) const { return calculateScalarImpl<double>(par); }

    void calculateAt(const ResultRequirementsType &req, const Eigen::Vector<double, Traits::mydim> &local,
                     ResultTypeMap<double> &result) const {
      using namespace Dune::Indices;
    }

   protected:
    template<typename ScalarType>
    auto calculateScalarImpl(const FERequirementType &par, const std::optional<const Eigen::VectorX<ScalarType>> &dx
    = std::nullopt) const -> ScalarType {
      using namespace Dune::Indices;
      ScalarType energy = 0.0;

            const auto &lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      using namespace Dune::DerivativeDirections;
      using namespace Dune::Indices;
      const auto [ displacementFunction, directorFunction,
                  directorReferenceFunction]
          = createFunctions(par,dx);

      KinematicVariables<ScalarType> kin{};

      Dune::BlockVector<Dune::FieldMatrix<ScalarType, 8, 3>> bopMidSurface(numNodeMidSurface);
      Dune::BlockVector<Dune::FieldMatrix<ScalarType, 8, 2>> bopDirector(numNodeDirector);
      const int midSurfaceDofs = numNodeMidSurface * midSurfaceDim;
      for (const auto &[gpIndex, gp] : displacementFunction.viewOverIntegrationPoints()) {
        const auto weight = geo_->integrationElement(gp.position()) * gp.weight();

        kin.t           = directorFunction.evaluate(gpIndex,Dune::on(referenceElement));
        kin.t0          = directorReferenceFunction.evaluate(gpIndex,Dune::on(referenceElement));
        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
        kin.A1andA2     = Dune::toEigen(geo_->jacobianTransposed(gp.position())).transpose();
        kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));
        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll),Dune::on(referenceElement));

        auto Egl                = calculateGreenLagrangianStrains(kin);
        const auto StimesWeight = (calculateStressResultants(CMat_, Egl)).eval();
        energy += 0.5 * Egl.dot(StimesWeight) * weight;
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : displacementFunction.viewOverIntegrationPoints()) {
          const auto uVal                              = displacementFunction.evaluate(gpIndex);
          Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(geo_->global(gp.position()), lambda)[0];
          energy -= uVal.dot(fext) * geo_->integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary) return energy;

      auto element = localViewFlat_.element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, 1>::rule(intersection.type(), displacementFunction.order());

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, 2>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto uVal = displacementFunction.evaluate(quadPos);

          // Value of the Neumann data at the current position
          auto neumannValue = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda)[0];

          energy -= neumannValue.dot(uVal) * curQuad.weight() * integrationElement;
        }
      }
      return energy;
    }

    mutable Eigen::Matrix<double, Eigen::Dynamic, Traits::mydim> dNA;

    const MultiTypeVector &x0;

    LocalViewBlocked localViewBlocked_;
    LocalViewFlat localViewFlat_;
    Dune::CachedLocalBasis<std::remove_cvref_t<
        decltype(std::declval<LocalViewFlat>().tree().child(Dune::Indices::_0, 0).finiteElement().localBasis())>>
        localBasisMidSurface;
    Dune::CachedLocalBasis<std::remove_cvref_t<
        decltype(std::declval<LocalViewFlat>().tree().child(Dune::Indices::_1, 0).finiteElement().localBasis())>>
        localBasisDirector;

    int order{};
    using LoadFunctorType =
    std::function<std::array<Eigen::Vector<double, Traits::worlddim>,2>(const Dune::FieldVector<double, Traits::worlddim>&,
                                                          const double&)>;
    LoadFunctorType  volumeLoad;
    LoadFunctorType  neumannBoundaryLoad;

    const BoundaryPatch<GridView>* neumannBoundary;
    std::shared_ptr<const Geometry> geo_;

    YoungsModulusAndPoissonsRatio material;
    double thickness_;
    Eigen::Matrix<double, 8, 8> CMat_;

    size_t numNodeMidSurface{};
    size_t numNodeDirector{};
  };

}  // namespace Ikarus

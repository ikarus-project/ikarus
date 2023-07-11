// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>
#include <iosfwd>
#include <numbers>

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
#include "membranestrains.hh"
#include <ikarus/utils/tensorUtils.hh>
#include "directorFunctions.hh"

namespace Ikarus {
template<typename ResultantBasedShell>
struct StressBasedShellRM;

  template <typename Basis_, bool useEigenRef_ = false>
  class NonLinearRMshell {
   public:
    using Basis = Basis_;

    static constexpr bool isEuclidean = false;

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

    using FERequirementType      = FErequirements<std::reference_wrapper<MultiTypeVector>>;
    using ResultRequirementsType = ResultRequirements<MultiTypeVector>;
    using LocalViewBlocked       = typename Basis::UntouchedBasis::LocalView;
    using LocalViewFlat          = typename Basis::FlatBasis::LocalView;
    using Geometry               = typename LocalViewFlat::Element::Geometry;
    using GridView          = typename Basis::FlatBasis::GridView;
    using Traits = TraitsFromLocalView<LocalViewFlat, false>;
    static constexpr int useEigenRef = useEigenRef_;

    static constexpr int myDim = Traits::mydim;
    static constexpr int worlddim = Traits::worlddim;
    mutable MembraneStrain membraneStrain;
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
           const auto &simulationFlag = fESettings.request<int>("simulationFlag");

           if (simulationFlag==0)
             membraneStrain = DefaultMembraneStrain();
           else if (simulationFlag==1)
             membraneStrain = CASMembraneStrain<CASAnsatzFunction>();
           else if (simulationFlag==2)
             membraneStrain = CASMembraneStrain<CASAnsatzFunctionANS>();

      CMat_.setZero();
//      // membrane
//      const double fac1 = thickness_ * Emodul / (1 - nu * nu);
//      CMat_(0, 0) = CMat_(1, 1) = fac1;
//      CMat_(2, 2)               = fac1 * (1 - nu) * 0.5;
//      CMat_(1, 0) = CMat_(0, 1) = fac1 * nu;
//
//      // bending
//      const double fac2 = thickness_ * thickness_ * thickness_ / 12 * Emodul / (1 - nu * nu);
//      CMat_(3, 3) = CMat_(4, 4) = fac2;
//      CMat_(5, 5)               = fac2 * (1 - nu) * 0.5;
//      CMat_(3, 4) = CMat_(4, 3) = fac2 * nu;
//
//      // transverse shear
//      const double fac3 = kappa_ * thickness_ * Emodul * 0.5 / (1 + nu);
//      CMat_(6, 6) = CMat_(7, 7) = fac3;

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

    template<typename ScalarType>
    auto calc3DMetric(const KinematicVariables<ScalarType> &kin, double zeta) const
    {
      const auto G1 = (kin.A1() + zeta*kin.t0d1()).eval();
      const auto G2 =(kin.A2() + zeta*kin.t0d2()).eval();
      const auto G3 = kin.t0;
      Eigen::Matrix<double,3,3> J3D;
      J3D.col(0)=G1;
      J3D.col(1)=G2;
      J3D.col(2)=G3;

      Eigen::Matrix<double,3,3> G= J3D.transpose()*J3D; //here quadratic parts of zeta appear!
      return G;
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
    auto computeMaterialAndStrains(const Dune::FieldVector<double, 2> &gpPos,
                                                const Geometry &geo,
                                                const auto &uFunction,const KinematicVariables<ScalarType> &kin) const {
      Eigen::Vector3<ScalarType> epsV=membraneStrain.value(gpPos,geo,uFunction);
      const Eigen::Matrix<double, 2, 2> A = kin.A1andA2.transpose()*kin.A1andA2;
      Eigen::Matrix<double, 3, 3> G;
      G.setZero();
      G.block<2, 2>(0, 0) = A;
      G(2, 2) = 1;
      const Eigen::Matrix<double, 3, 3> GInv = G.inverse();
      const auto C = materialTangent(GInv);
      // membrane
      Eigen::Vector3<ScalarType> kappaV;

      kappaV<<  kin.a1().dot(kin.td1()) - kin.A1().dot(kin.t0d1()), kin.a2().dot(kin.td2()) - kin.A2().dot(kin.t0d2()),
          kin.a1().dot(kin.td2()) + kin.a2().dot(kin.td1()) - kin.A1().dot(kin.t0d2()) - kin.A2().dot(kin.t0d1());
      Eigen::Vector2<ScalarType> gammaV;

      gammaV<<kin.a1().dot(kin.t) - kin.A1().dot(kin.t0), kin.a2().dot(kin.t) - kin.A2().dot(kin.t0);

      return std::make_tuple(C,epsV,kappaV,gammaV);
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

      for (int i = 0; i < nDofs1; i++) {
        int localIndexI = 0;
//        if (i < nDofs0) {
//          auto &node  = localViewBlocked_.tree().child(_0, 0);
//          localIndexI = node.localIndex(i);
//        } else {
          auto &node  = localViewBlocked_.tree().child(_1);
          localIndexI = node.localIndex(i );
//        }
        auto multiIndex = localViewBlocked_.index(localIndexI);

        // The CompositeBasis number is contained in multiIndex[0]
        // multiIndex[1] contains the actual index
//        if (multiIndex[0] == 0)
//          localConfiguration0[i] = x0[_0][multiIndex[1]];
//        else
          if (multiIndex[0] == 1)
          localConfiguration1[i] = x0[_1][multiIndex[1]];
      }
      return localConfiguration1;
    }

//    template<typename ScalarType>
//    Eigen::Matrix<ScalarType, 8, 3> boperatorMidSurface(const KinematicVariables<ScalarType> &kin, int integrationPointIndex,
//                                                        int coeffIndex, const auto &displacementFunction) const {
//      using namespace Dune::TypeTree::Indices;
//      using namespace Dune::DerivativeDirections;
//      const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda2
//          = displacementFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll, coeff(_0, coeffIndex)),Dune::on(referenceElement));
//      const auto &diffa1 = diffa1Anda2[0];
//      const auto &diffa2 = diffa1Anda2[1];
//      Eigen::Matrix<ScalarType, 8, 3> bop;
//      bop.setZero();
//      bop.row(0) = kin.a1().transpose() * diffa1;  // membrane_{,disp}
//      bop.row(1) = kin.a2().transpose() * diffa2;
//      bop.row(2) = kin.a2().transpose() * diffa1 + kin.a1().transpose() * diffa2;
//
//      bop.row(3) = kin.td1().transpose() * diffa1;  // bending_{,disp}
//      bop.row(4) = kin.td2().transpose() * diffa2;
//      bop.row(5) = kin.td2().transpose() * diffa1 + kin.td1().transpose() * diffa2;
//
//      bop.row(6) = kin.t.transpose() * diffa1;  // trans_shear_{,disp}
//      bop.row(7) = kin.t.transpose() * diffa2;
//      return bop;
//    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 3, 3> boperatorMidSurfaceBending(const KinematicVariables<ScalarType> &kin, int integrationPointIndex,
                                                        int coeffIndex, const auto &displacementFunction) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda2
          = displacementFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll, coeff(_0, coeffIndex)),Dune::on(referenceElement));
      const auto &diffa1 = diffa1Anda2[0];
      const auto &diffa2 = diffa1Anda2[1];
      Eigen::Matrix<ScalarType, 3, 3> bop;
      bop.setZero();
      bop.row(0) = kin.td1().transpose() * diffa1;  // bending_{,disp}
      bop.row(1) = kin.td2().transpose() * diffa2;
      bop.row(2) = kin.td2().transpose() * diffa1 + kin.td1().transpose() * diffa2;

      return bop;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 3, 2> kgMidSurfaceDirectorBending(const KinematicVariables<ScalarType> &kin,const auto& N,const auto& dN, int integrationPointIndex,
                                                               int i,int j, const auto &displacementFunction, const auto &directorFunction,const auto& S) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const double& Ni   = N[i];
      const double& dN1i = dN(i, 0);
      const double& dN2i = dN(i, 1);
      const std::array<Eigen::Matrix<ScalarType, 3, 2>, 2> WJ
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll, coeff(_1, j)),Dune::on(referenceElement));
      Eigen::Matrix<ScalarType, 3, 2> kg = S[3] * dN1i * WJ[0] + S[4] * dN2i * WJ[1] + S[5] * (dN1i * WJ[1] + dN2i * WJ[0]);  // bending_{,dir,disp}*M

      return kg;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 2, 2> kgDirectorDirectorBending(const KinematicVariables<ScalarType> &kin,const auto& Nd,const auto& dNd, int integrationPointIndex,
                                                                int i,int j, const auto &displacementFunction, const auto &directorFunction,const auto& S) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const double& Ni   = Nd[i];
      const double& dN1i = dNd(i, 0);
      const double& dN2i = dNd(i, 1);
      const double& Nj = Nd[j];
      const double& dN1j = dNd(j, 0);
      const double& dN2j = dNd(j, 1);

      const double NdN1 = dN1j * Ni + Nj * dN1i;
      const double NdN2 = dN2j * Ni + Nj * dN2i;
//      const Eigen::Vector<ScalarType,3> SM = S.template segment<3>(3);
      const auto a1 = kin.a1().eval();
      const auto a2 = kin.a2().eval();
      const Eigen::Matrix<ScalarType, 2, 2> fac11
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatial(0),coeff(_1, i,_1, j)),Dune::along(a1),Dune::on(referenceElement));
      const Eigen::Matrix<ScalarType, 2, 2> fac21
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatial(1),coeff(_1, i,_1, j)),Dune::along(a1),Dune::on(referenceElement));
      const Eigen::Matrix<ScalarType, 2, 2> fac12
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatial(0),coeff(_1, i,_1, j)),Dune::along(a2),Dune::on(referenceElement));
      const Eigen::Matrix<ScalarType, 2, 2> fac22
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatial(1),coeff(_1, i,_1, j)),Dune::along(a2),Dune::on(referenceElement));
      Eigen::Matrix<ScalarType, 2, 2> kg = fac11 * S[3] + fac22 * S[4] +(fac12 + fac21) * S[5];  // bending_{,dir,dir}*M

      return kg;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 2, 2> kgSecondDirectorDirectorBending(const KinematicVariables<ScalarType> &kin,const auto& Nd,const auto& dNd, int integrationPointIndex,
                                                              int i,int j, const auto &displacementFunction, const auto &directorFunction,const auto& S) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const double& Ni   = Nd[i];
      const double& dN1i = dNd(i, 0);
      const double& dN2i = dNd(i, 1);
      const double& Nj = Nd[j];
      const double& dN1j = dNd(j, 0);
      const double& dN2j = dNd(j, 1);

      const double NdN1 = dN1j * Ni + Nj * dN1i;
      const double NdN2 = dN2j * Ni + Nj * dN2i;
//      const Eigen::Vector<ScalarType,3> SM = S.template segment<3>(3);
      const auto td1 = kin.td1().eval();
      const auto td2 = kin.td2().eval();
      const Eigen::Matrix<ScalarType, 2, 2> fac11
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatial(0),coeff(_1, i,_1, j)),Dune::along(td1),Dune::on(referenceElement));
      const Eigen::Matrix<ScalarType, 2, 2> fac21
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatial(1),coeff(_1, i,_1, j)),Dune::along(td1),Dune::on(referenceElement));
      const Eigen::Matrix<ScalarType, 2, 2> fac12
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatial(0),coeff(_1, i,_1, j)),Dune::along(td2),Dune::on(referenceElement));
      const Eigen::Matrix<ScalarType, 2, 2> fac22
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatial(1),coeff(_1, i,_1, j)),Dune::along(td2),Dune::on(referenceElement));

      const auto WI
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll,coeff(_1, i)),Dune::on(referenceElement));
      const auto WJ
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll,coeff(_1, j)),Dune::on(referenceElement));
      Eigen::Matrix<ScalarType, 2, 2> kg;
      kg  = WI[0].transpose() * WJ[0] * S[0] + WI[1].transpose() * WJ[1] * S[1] + (WI[0] .transpose()* WJ[1] + WI[1].transpose() * WJ[0]) * S[2]; //{/zeta^2* t_{,a}\cdot t_{,a}}_{,dir,dir}
      kg += fac11 * S[0] + fac22 * S[1] + (fac21 + fac12) * S[2];

      return kg;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 2, 2> kgDirectorDirectorShear(const KinematicVariables<ScalarType> &kin,const auto& Nd,const auto& dNd, int integrationPointIndex,
                                                                 int i,int j, const auto &displacementFunction, const auto &directorFunction,const auto& S) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const double& Ni   = Nd[i];
      const double& dN1i = dNd(i, 0);
      const double& dN2i = dNd(i, 1);
      const double& Nj = Nd[j];
      const double& dN1j = dNd(j, 0);
      const double& dN2j = dNd(j, 1);

      const double NdN1 = dN1j * Ni + Nj * dN1i;
      const double NdN2 = dN2j * Ni + Nj * dN2i;
//      const Eigen::Vector<ScalarType,3> SM = S.template segment<3>(3);
      const auto a1 = kin.a1().eval();
      const auto a2 = kin.a2().eval();
      const Eigen::Matrix<ScalarType, 2, 2> S1
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(coeff(_1, i,_1, j)),Dune::along(a1),Dune::on(referenceElement));

      const Eigen::Matrix<ScalarType, 2, 2> S2
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(coeff(_1, i,_1, j)),Dune::along(a2),Dune::on(referenceElement));
      Eigen::Matrix<ScalarType, 2, 2> kg = S1* S[6] + S2 * S[7];

      return kg;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 3, 2> kgMidSurfaceDirectorShear(const KinematicVariables<ScalarType> &kin,const auto& N,const auto& dN, int integrationPointIndex,
                                                                int i,int j, const auto &displacementFunction, const auto &directorFunction,const auto& S) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const double& Ni   = N[i];
      const double& dN1i = dN(i, 0);
      const double& dN2i = dN(i, 1);
      const Eigen::Matrix<ScalarType, 3, 2> P
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(coeff(_1, j)),Dune::on(referenceElement));
      Eigen::Matrix<ScalarType, 3, 2> kg = P * (dN1i * S[6] + dN2i * S[7]); // shear_{,dir,disp}*Q
      // bending_{,dir,disp}*M

      return kg;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 2, 3> boperatorMidSurfaceShear(const KinematicVariables<ScalarType> &kin, int integrationPointIndex,
                                                        int coeffIndex, const auto &displacementFunction) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda2
          = displacementFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll, coeff(_0, coeffIndex)),Dune::on(referenceElement));
      const auto &diffa1 = diffa1Anda2[0];
      const auto &diffa2 = diffa1Anda2[1];
      Eigen::Matrix<ScalarType, 2, 3> bop;
      bop.setZero();
      bop.row(0) = kin.t.transpose() * diffa1;  // trans_shear_{,disp}
      bop.row(1) = kin.t.transpose() * diffa2;
      return bop;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 3, 2> boperatorDirectorBending(const Eigen::Matrix<ScalarType, 3, 2>& j, int integrationPointIndex,
                                                      int coeffIndex, const auto &directorFunction) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const std::array<Eigen::Matrix<ScalarType, 3, 2>, 2> diffdt
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll, coeff(_1, coeffIndex)),Dune::on(referenceElement));
      const Eigen::Matrix<ScalarType, 3, 2> difft
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(coeff(_1, coeffIndex)),Dune::on(referenceElement));

      const auto &diffdtd1 = diffdt[0];
      const auto &diffdtd2 = diffdt[1];

      Eigen::Matrix<ScalarType, 3, 2> bop;
      bop.setZero();
      bop.row(0) = j.col(0).transpose() * diffdtd1;  // bending_{,disp}
      bop.row(1) = j.col(1).transpose() * diffdtd2;
      bop.row(2) = j.col(1).transpose() * diffdtd1 + j.col(0).transpose() * diffdtd2;

      return bop;
    }

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, 2, 2> boperatorDirectorShear(const KinematicVariables<ScalarType> &kin, int integrationPointIndex,
                                                      int coeffIndex, const auto &directorFunction) const {
      using namespace Dune::TypeTree::Indices;
      using namespace Dune::DerivativeDirections;
      const std::array<Eigen::Matrix<ScalarType, 3, 2>, 2> diffdt
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(spatialAll, coeff(_1, coeffIndex)),Dune::on(referenceElement));
      const Eigen::Matrix<ScalarType, 3, 2> difft
          = directorFunction.evaluateDerivative(integrationPointIndex, Dune::wrt(coeff(_1, coeffIndex)),Dune::on(referenceElement));

      const auto &diffdtd1 = diffdt[0];
      const auto &diffdtd2 = diffdt[1];

      Eigen::Matrix<ScalarType, 2, 2> bop;

      bop.row(0) = kin.a1().transpose() * difft;  // trans_shear_{,disp}
      bop.row(1) = kin.a2().transpose() * difft;

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
          const auto globalIndex = localViewBlocked_.index(localViewBlocked_.tree().child(_0, 0).localIndex(i));
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            displacements[i][k2] = mNodal[globalIndex[1]][k2] +dx.value()[i*worlddim + k2];
        }
        const int nDofs0 = this->localViewBlocked().tree().child(_0, 0).finiteElement().size()*3;

        for (auto i = 0U; i < fe1.size(); ++i) {
          const int localIndex = localViewBlocked_.tree().child(_1, 0).localIndex(i);
          const auto globalIndex = localViewBlocked_.index(localIndex);
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            localDirectorConfiguration[i][k2] = dNodal[globalIndex[1]][k2]+dx.value()[nDofs0+i*worlddim + k2];
        }
      }else
      {        for (auto i = 0U; i < fe0.size(); ++i) {
          const auto globalIndex = localViewBlocked_.index(localViewBlocked_.tree().child(_0).localIndex(i));
          displacements[i] = mNodal[globalIndex[1]];
        }

        for (auto i = 0U; i < fe1.size(); ++i) {
          const int localIndex = localViewBlocked_.tree().child(_1, 0).localIndex(i);
          const auto globalIndex = localViewBlocked_.index(localIndex);
          localDirectorConfiguration[i] = dNodal[globalIndex[1]];
        }
      }

      const auto localRefDirectorConfiguration = getReferenceLocalConfigurations();


      using namespace Dune::DerivativeDirections;

      Dune::StandardLocalFunction displacementFunction(localBasisMidSurface, displacements, geo_, _0);
      Dune::ProjectionBasedLocalFunction2 directorFunctionImpl(localBasisDirector, localDirectorConfiguration, geo_, _1);
      Dune::ProjectionBasedLocalFunction2 directorReferenceFunctionImpl(localBasisDirector, localRefDirectorConfiguration,
                                                                   geo_, _1);
      using DirectorCurType= decltype(directorFunctionImpl);
      using DuneBasis = typename  DirectorCurType::DuneBasis;
      using CoeffContainer = typename  DirectorCurType::CoeffContainer;
      using GeometryL = typename  DirectorCurType::Geometry;
      static constexpr int orderID = DirectorCurType::id[0];
      using DirVariantCur = DirectorFunctionVar<DuneBasis, CoeffContainer,GeometryL,orderID>;
      DirVariantCur directorFunction(directorFunctionImpl);

      using DirectorRefType= decltype(directorReferenceFunctionImpl);
      using DuneBasisRef = typename  DirectorRefType::DuneBasis;
      using CoeffContainerRef = typename  DirectorRefType::CoeffContainer;
      using GeometryLRef = typename  DirectorRefType::Geometry;
      static constexpr int orderIDRef = DirectorRefType::id[0];
      using DirVariantRef = DirectorFunctionVar<DuneBasisRef, CoeffContainerRef,GeometryLRef,orderIDRef>;
      DirVariantRef directorReferenceFunction(directorReferenceFunctionImpl);
      return std::make_tuple( displacementFunction, directorFunction,
                             directorReferenceFunction);
    }

    inline void calculateMatrix(const FERequirementType &par, typename Traits::template MatrixType<> K) const {
      calculateMatrixImpl<double>(par, K);
    }
    template<typename ScalarType>
    void calculateMatrixImpl(const FERequirementType &par, typename Traits::template MatrixType<ScalarType> hred,
                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
      using namespace Dune::Indices;
      hred.setZero();

      //      const auto &lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      using namespace Dune::DerivativeDirections;
      using namespace Dune::Indices;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = createFunctions(par,dx);

      KinematicVariables<ScalarType> kin{};

      Eigen::Matrix<ScalarType, 8, 3> bopMidSurfaceI;
      Eigen::Matrix<ScalarType, 8, 2> bopDirectorI;

      Eigen::Matrix<ScalarType, 8, 3> bopMidSurfaceJ;
      Eigen::Matrix<ScalarType, 8, 2> bopDirectorJ;

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

        auto [C3D,epsV,kappV,gammaV]                = computeMaterialAndStrains(gp.position(),*geo_,displacementFunction,kin);
        Eigen::Vector<ScalarType,8> Egl,S;
        Egl<<epsV,kappV,gammaV;
        S.template segment<3>(0)= thickness_*C3D.template block<3,3>(0,0)*epsV;
        S.template segment<3>(3)= Dune::power(thickness_, 3)/12.0*C3D.template block<3,3>(0,0)*kappV;
        S.template segment<2>(6)= thickness_*C3D.template block<2,2>(3,3)*gammaV;
        const Eigen::Matrix<ScalarType,2,3> jE = kin.a1anda2.transpose();
        CMat_. template block<3,3>(0,0)=thickness_*C3D.template block<3,3>(0,0);
        CMat_. template block<3,3>(3,3)=Dune::power(thickness_, 3)/12.0*C3D.template block<3,3>(0,0);
        CMat_. template block<2,2>(6,6)=thickness_*C3D.template block<2,2>(3,3);
        const auto &Nd = localBasisMidSurface.evaluateJacobian(gpIndex);
        const auto &N = localBasisMidSurface.evaluateFunction(gpIndex);
        const auto &dNdirector = localBasisDirector.evaluateJacobian(gpIndex);
        const auto &Ndirector = localBasisDirector.evaluateFunction(gpIndex);

        for (int i = 0; i < numNodeMidSurface; ++i) {
          const auto indexI = midSurfaceDim * i;
          const auto bopIMembraneI   = membraneStrain.derivative(gp.position(),jE, Nd, *geo_,displacementFunction,localBasisMidSurface, i);
          const auto bopIBendingI   = boperatorMidSurfaceBending(kin,gpIndex, i,displacementFunction);
          const auto bopIShearI   = boperatorMidSurfaceShear(kin,gpIndex, i,displacementFunction);
          bopMidSurfaceI<< bopIMembraneI,bopIBendingI, bopIShearI;
          for (int j = i; j < numNodeMidSurface; ++j) {
            const auto indexJ = midSurfaceDim * j;

            const auto bopIMembraneJ   = membraneStrain.derivative(gp.position(),jE, Nd, *geo_,displacementFunction,localBasisMidSurface, j);
            const auto bopIBendingJ   = boperatorMidSurfaceBending(kin,gpIndex, j,displacementFunction);
            const auto bopIShearJ   = boperatorMidSurfaceShear(kin,gpIndex, j,displacementFunction);
            bopMidSurfaceJ<< bopIMembraneJ,bopIBendingJ, bopIShearJ;

            hred.template block<midSurfaceDim, midSurfaceDim>(indexI, indexJ)
                += bopMidSurfaceI.transpose() * CMat_ * bopMidSurfaceJ * weight;

            Eigen::Matrix<ScalarType, 3, 3> kgMembraneIJ = membraneStrain.secondDerivative(gp.position(),Nd,*geo_,displacementFunction,localBasisMidSurface, S.template segment<3>(0).eval(), i, j);
            hred.template block<3, 3>(3*i, 3*j) += kgMembraneIJ*weight;
          }

          for (int j = 0; j < numNodeDirector; ++j) {
            const auto indexJ = midSurfaceDofs + directorCorrectionDim * j;
            const auto bopBendingJ   = boperatorDirectorBending(kin.a1anda2, gpIndex, j, directorFunction);
            const auto bopShearJ   = boperatorDirectorShear(kin, gpIndex, j, directorFunction);

            bopDirectorJ <<Eigen::Matrix<double,3,2>::Zero(),bopBendingJ,bopShearJ;

            Eigen::Matrix<ScalarType, 3, 2> kg= kgMidSurfaceDirectorBending(kin,N,Nd,gpIndex,i,j,displacementFunction,directorFunction,S);
            Eigen::Matrix<ScalarType, 3, 2> kg2= kgMidSurfaceDirectorShear(kin,N,Nd,gpIndex,i,j,displacementFunction,directorFunction,S);

            hred.template block<midSurfaceDim, directorCorrectionDim>(indexI, indexJ)
                += bopMidSurfaceI.transpose() * CMat_ * bopDirectorJ * weight;

            hred.template block<midSurfaceDim, directorCorrectionDim>(indexI, indexJ)
                += (kg+kg2) * weight;
          }
        }

        for (int i = 0; i < numNodeDirector; ++i) {
          const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
          const auto bopBendingI   = boperatorDirectorBending(kin.a1anda2, gpIndex, i, directorFunction);
          const auto bopShearI   = boperatorDirectorShear(kin, gpIndex, i, directorFunction);

          bopDirectorI <<Eigen::Matrix<double,3,2>::Zero(),bopBendingI,bopShearI;

          for (int j = i; j < numNodeDirector; ++j) {
            const auto indexJ = midSurfaceDofs + directorCorrectionDim * j;

            const auto bopBendingJ   = boperatorDirectorBending(kin.a1anda2, gpIndex, j, directorFunction);
            const auto bopShearJ   = boperatorDirectorShear(kin, gpIndex, j, directorFunction);
            bopDirectorJ <<Eigen::Matrix<double,3,2>::Zero(),bopBendingJ,bopShearJ;

            Eigen::Matrix<ScalarType, 2, 2> kgBending= kgDirectorDirectorBending(kin,Ndirector,dNdirector,gpIndex,i,j,displacementFunction,directorFunction,S);
            Eigen::Matrix<ScalarType, 2, 2> kgShear= kgDirectorDirectorShear(kin,Ndirector,dNdirector,gpIndex,i,j,displacementFunction,directorFunction,S);


            hred.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ)
                += bopDirectorI.transpose() * CMat_ * bopDirectorJ * weight;

            hred.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ)
                += (kgBending+kgShear)* weight;
          }
        }
      }
      hred.template triangularView<Eigen::StrictlyLower>() = hred.transpose();
    }

    //    [[nodiscard]] int size() const { return localView.size(); }

    inline void calculateVector(const FERequirementType &par, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(par, force);
    }
    template<typename ScalarType>
    void calculateVectorImpl(const FERequirementType &par, typename Traits::template VectorType<ScalarType> rieGrad,
                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
      using namespace Dune::Indices;
      rieGrad.setZero();

      const auto &lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      using namespace Dune::DerivativeDirections;
      using namespace Dune::Indices;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = createFunctions(par,dx);

      KinematicVariables<ScalarType> kin{};

      Eigen::Matrix<ScalarType, 8, 3> bopMidSurface;
      Eigen::Matrix<ScalarType, 8, 2> bopDirector;
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

        auto [C3D,epsV,kappV,gammaV]                = computeMaterialAndStrains(gp.position(),*geo_,displacementFunction,kin);
        Eigen::Vector<ScalarType,8> Egl,S;
        Egl<<epsV,kappV,gammaV;
        S.template segment<3>(0)= thickness_*C3D.template block<3,3>(0,0)*epsV;
        S.template segment<3>(3)= Dune::power(thickness_, 3)/12.0*C3D.template block<3,3>(0,0)*kappV;
        S.template segment<2>(6)= thickness_*C3D.template block<2,2>(3,3)*gammaV;
        const Eigen::Matrix<ScalarType,2,3> jE = kin.a1anda2.transpose();
        const auto &Nd = localBasisMidSurface.evaluateJacobian(gpIndex);
        for (int i = 0; i < numNodeMidSurface; ++i) {
          const auto indexI = midSurfaceDim * i;
          const auto bopIMembrane   = membraneStrain.derivative(gp.position(),jE, Nd, *geo_,displacementFunction,localBasisMidSurface, i);
          const auto bopIBending   = boperatorMidSurfaceBending(kin,gpIndex, i,displacementFunction);
          const auto bopIShear   = boperatorMidSurfaceShear(kin,gpIndex, i,displacementFunction);
          bopMidSurface<< bopIMembrane,bopIBending, bopIShear;
          rieGrad.template segment<midSurfaceDim>(indexI) += bopMidSurface.transpose() * S*weight;
        }

        for (int i = 0; i < numNodeDirector; ++i) {
          const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
          const auto bopIBending   = boperatorDirectorBending(kin.a1anda2, gpIndex, i, directorFunction);
          const auto bopIShear   = boperatorDirectorShear(kin, gpIndex, i, directorFunction);

          bopDirector<<Eigen::Matrix<double,3,2>::Zero(),bopIBending,bopIShear;
          rieGrad.template segment<directorCorrectionDim>(indexI) += bopDirector.transpose() * S*weight;
        }
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : displacementFunction.viewOverIntegrationPoints()) {
          const auto [fext,mext] = volumeLoad(geo_->global(gp.position()), lambda);
          for (size_t i = 0; i < numNodeMidSurface; ++i) {
            const auto indexI = midSurfaceDim * i;
            const auto udCi = displacementFunction.evaluateDerivative(gpIndex, Dune::wrt(coeff(_0,i)));

            rieGrad.template segment<midSurfaceDim>(indexI)
                -= (udCi * fext) * geo_->integrationElement(gp.position()) * gp.weight();
          }
          for (size_t i = 0; i < numNodeDirector; ++i) {
            const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
            const auto tdCi = directorFunction.evaluateDerivative(gpIndex, Dune::wrt(coeff(_1,i)));
            rieGrad.template segment<directorCorrectionDim>(indexI)
                -= (tdCi.transpose() * mext) * geo_->integrationElement(gp.position()) * gp.weight();
          }
        }
      }

      // External forces, boundary forces, i.e. at the Neumann boundary
      if (not neumannBoundary) return;

      auto element = localViewFlat_.element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary->contains(intersection)) continue;

        // Integration rule along the boundary
        const auto& quadLine = Dune::QuadratureRules<double, 1>::rule(intersection.type(), displacementFunction.order());

        for (const auto& curQuad : quadLine) {
          const Dune::FieldVector<double, 2>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());
          const auto [fext,mext]
              = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda);
          // The value of the local function wrt the i-th coef
          for (size_t i = 0; i < numNodeMidSurface; ++i) {
            const auto indexI = midSurfaceDim * i;
            const auto udCi = displacementFunction.evaluateDerivative(quadPos, Dune::wrt(coeff(_0,i)));

            rieGrad.template segment<midSurfaceDim>(indexI) -= udCi * fext * curQuad.weight() * integrationElement;
          }

          for (size_t i = 0; i < numNodeDirector; ++i) {
            const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
            const auto tdCi = directorFunction.evaluateDerivative(quadPos, Dune::wrt(coeff(_1,i)));
            rieGrad.template segment<directorCorrectionDim>(indexI)
                -= (tdCi.transpose() * mext) * curQuad.weight() * integrationElement;;
          }
        }
      }
    }

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

        auto [C3D,epsV,kappV,gammaV]                = computeMaterialAndStrains(gp.position(),*geo_,displacementFunction,kin);
        Eigen::Vector<ScalarType,8> Egl,S;
        Egl<<epsV,kappV,gammaV;
        S.template segment<3>(0)= thickness_*C3D.template block<3,3>(0,0)*epsV;
        S.template segment<3>(3)= Dune::power(thickness_, 3)/12.0*C3D.template block<3,3>(0,0)*kappV;
        S.template segment<2>(6)= thickness_*C3D.template block<2,2>(3,3)*gammaV;
        energy += 0.5 * Egl.dot(S) * weight;
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

    auto materialTangent(const Eigen::Matrix<double, 3, 3> &Aconv) const {
      const auto &emod_ = fESettings.request<double>("youngs_modulus");
      const auto &nu_ = fESettings.request<double>("poissons_ratio");
      const double lambda = emod_*nu_/((1.0 + nu_)*(1.0 - 2.0*nu_));
      const double mu = emod_/(2.0*(1.0 + nu_));
      const double lambdbar = 2.0*lambda*mu/(lambda + 2.0*mu);
      Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>> moduli;
      const auto AconvT = TensorCast(Aconv, std::array<Eigen::Index, 2>({3, 3}));
      moduli = lambdbar*dyadic(AconvT, AconvT).eval() + 2.0*mu*symmetricFourthOrder<double>(Aconv, Aconv);

      auto C = toVoigt(moduli);
      Eigen::Matrix<double, 5, 5> midSurf = C({0, 1, 5,4, 3}, {0, 1, 5,4, 3});
      return midSurf;
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
    mutable Eigen::Matrix<double, 8, 8> CMat_;
    friend StressBasedShellRM<NonLinearRMshell>;
    size_t numNodeMidSurface{};
    size_t numNodeDirector{};

    const LocalViewFlat& localView() const { return localViewFlat_; }
    LocalViewFlat& localView() { return localViewFlat_; }
  };

}  // namespace Ikarus

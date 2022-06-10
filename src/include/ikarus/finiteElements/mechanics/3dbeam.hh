
// /*
//  *  This file is part of the Ikarus distribution (https://github.com/rath3t/Ikarus).
//  *  Copyright (c) 2021 Alexander Müller.
//  *  Institut fuer Baustatik und Baudynamik
//  *  Universität Stuttgart
//  *
//  *  This library is free software; you can redistribute it and/or
//  *   modify it under the terms of the GNU Lesser General Public
//  *   License as published by the Free Software Foundation; either
//  *   version 2.1 of the License, or (at your option) any later version.
//
// *   This library is distributed in the hope that it will be useful,
// *   but WITHOUT ANY WARRANTY; without even the implied warranty of
// *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// *   Lesser General Public License for more details.
//
// *   You should have received a copy of the GNU Lesser General Public
// *   License along with this library; if not, write to the Free Software
// *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// *  USA
// *

#pragma once
#include "src/include/ikarus/finiteElements/feTraits.hh"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <concepts>
#include <iosfwd>
#include <numbers>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/localFunctions/impl/projectionBasedLocalFunction.hh>
#include <ikarus/localFunctions/impl/standardLocalFunction.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/manifolds/unitVector.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  namespace Impl {
    template <typename Derived>
    auto vecToSkewMatrix(const Eigen::MatrixBase<Derived>& vec) {
      constexpr int FieldSize = Derived::ColsAtCompileTime;
      Eigen::Matrix<typename Derived::Scalar, 3, FieldSize == 2 ? 1 : 3> skew;
      static_assert(Derived::RowsAtCompileTime == 1);

      if constexpr (FieldSize == 2) {
        skew[0] = vec(1);
        skew[1] = -vec(0);
        skew[2] = 0;
      } else if constexpr (FieldSize == 3) {
        skew << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
      }
      return skew;
    }

  }  // namespace Impl

  struct BeamMaterial {
    double E{0.0};
    double nu{0.0};
    double A{0.0};
    double I1{0.0};
    double I2{0.0};
    double J{0.0};
  };

  template<typename ScalarType>
  auto rotationMatrixColumnDerivatives(const Eigen::Vector<ScalarType,4>& q)
  {
    std::array<Eigen::Matrix<ScalarType,3,4> ,3> D;

    D[0] << q[0], -q[1], - q[2], q[3],
            q[1],  q[0],   q[3], q[2],
            q[2], -q[3],   q[0], - q[1];
    D[0]*=2.0;

    D[1] << q[1], q[0], - q[3], -q[2],
        -q[0],  q[1],   -q[2], q[3],
        q[3], q[2],   q[1], q[0];
    D[1]*=2.0;

    D[2] << q[2], q[3],  q[0], q[1],
        -q[3],  q[2],   q[1], -q[0],
        -q[0], -q[1],   q[2], q[3];
    D[2]*=2.0;

    return D;
  }

  template <typename BasisEmbedded, typename BasisReduced>
  class SimoReissnerBeam {
  public:
    static constexpr int centerLineDim           = BasisEmbedded::PreBasis::template SubPreBasis<0>::Node::CHILDREN;
    static constexpr int quaternionDim          = BasisEmbedded::PreBasis::template SubPreBasis<1>::Node::CHILDREN;
    static constexpr int quaternionCorrectionDim = quaternionDim - 1;
    using UnitQuaternionVector                       = Dune::BlockVector<Ikarus::UnitVector<double, quaternionDim>>;
    using CenterLinePositionsVector                      = Dune::BlockVector<Ikarus::RealTuple<double, centerLineDim>>;
    using MultiTypeVector                      = Dune::MultiTypeBlockVector<CenterLinePositionsVector, UnitQuaternionVector>;

    using FERequirementType      = FErequirements<MultiTypeVector>;
    using ResultRequirementsType = ResultRequirements<MultiTypeVector>;
    using LocalViewEmbedded      = typename BasisEmbedded::LocalView;
    using LocalViewReduced       = typename BasisReduced::LocalView;

    template <typename VolumeLoad>
    SimoReissnerBeam(BasisEmbedded& globalBasis, BasisReduced& globalBasisRed,
                                      const typename LocalViewEmbedded::Element& element, BeamMaterial& p_material,
                                      const VolumeLoad& p_volumeLoad)
        : localView_{globalBasis.localView()},
          localViewReduced{globalBasisRed.localView()},
          volumeLoad(p_volumeLoad),
          material{p_material}
    {
      localView_.bind(element);
      localViewReduced.bind(element);
      order         = 2 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();
      localBasisCenterLine
          = Ikarus::LocalBasis(localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis());
      localBasisCenterLine.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),
                         bindDerivatives(0, 1));
      localBasisQuaternions
          = Ikarus::LocalBasis(localView_.tree().child(Dune::Indices::_1, 0).finiteElement().localBasis());
      localBasisQuaternions.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),
                            bindDerivatives(0, 1));
    }

    using Traits = TraitsFromLocalView<LocalViewEmbedded>;

    using GlobalIndex = typename LocalViewReduced::MultiIndex;
    void globalIndices(std::vector<GlobalIndex>& globalIndices) const {
      using namespace Dune::Indices;
      const auto& fe = localViewReduced.tree().child(_0, 0).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < centerLineDim; ++j) {
          globalIndices.push_back(localViewReduced.index((localViewReduced.tree().child(_0, j).localIndex(i))));
        }
      }
      const auto& fe2 = localViewReduced.tree().child(_1, 0).finiteElement();
      for (size_t i = 0; i < fe2.size(); ++i) {
        for (int j = 0; j < quaternionCorrectionDim; ++j) {
          globalIndices.push_back(localViewReduced.index((localViewReduced.tree().child(_1, j).localIndex(i))));
        }
      }
    }

    mutable Eigen::MatrixXd hEuk;
    mutable Eigen::VectorXd eukGrad;
    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& hred) const {
      using namespace Dune::Indices;
      hEuk.resize(localView_.size(), localView_.size());
      eukGrad.resize(localView_.size());
      Eigen::VectorXdual2nd dx2nd;
      dx2nd.resize(localView_.size());
      dx2nd.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      autodiff::dual2nd e;
      autodiff::hessian(f, autodiff::wrt(dx2nd), at(dx2nd), e, eukGrad, hEuk);
      const auto& quaternions     = par.getSolution(Ikarus::FESolutions::displacementAndQuaternions)[_1];
      const auto& feDisp    = localView_.tree().child(_0, 0).finiteElement();
      const auto& feQuats = localView_.tree().child(_1, 0).finiteElement();

      const int dispHessianSize     = feDisp.size() * centerLineDim;
//      const int magFullHessianSize = feMag.size() * directorDim;

      for (auto i = 0U; i < feQuats.size(); ++i) {
        const size_t indexRedI  = dispHessianSize + i * quaternionCorrectionDim;
        const size_t indexI     = dispHessianSize + i * quaternionDim;
        const auto globalIndexI = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        const Eigen::Matrix<double, quaternionCorrectionDim, quaternionDim> BLAIT
            = quaternions[globalIndexI[1]].orthonormalFrame().transpose();
        for (auto j = 0U; j < feQuats.size(); ++j) {
          const size_t indexRedJ  = dispHessianSize + j * quaternionCorrectionDim;
          const size_t indexJ     = dispHessianSize + j * quaternionDim;
          const auto globalIndexJ = localView_.index(localView_.tree().child(_1, 0).localIndex(j));
          const auto BLAJ         = quaternions[globalIndexJ[1]].orthonormalFrame();

          hred.template block<quaternionCorrectionDim, quaternionCorrectionDim>(indexRedI, indexRedJ)
              = BLAIT * hEuk.block<quaternionDim, quaternionDim>(indexI, indexJ) * BLAJ;
        }
        hred.template block<quaternionCorrectionDim, quaternionCorrectionDim>(indexRedI, indexRedI)
            -= quaternions[globalIndexI[1]].getValue().dot(eukGrad.template segment<quaternionDim>(indexI))
               * Eigen::Matrix<double, quaternionCorrectionDim, quaternionCorrectionDim>::Identity();
      }



      for (auto i = 0U; i < feQuats.size(); ++i) {
        const size_t indexRedI  = dispHessianSize + i * quaternionCorrectionDim;
        const size_t indexI     = dispHessianSize + i * quaternionDim;
        const auto globalIndexI = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        const Eigen::Matrix<double, quaternionCorrectionDim, quaternionDim> BLAIT
            = quaternions[globalIndexI[1]].orthonormalFrame().transpose();
        for (auto j = 0U; j < feDisp.size(); ++j) {
          const size_t indexRedJ =  j * centerLineDim;
          const size_t indexJ    =  j * centerLineDim;
          hred.template block<quaternionCorrectionDim, centerLineDim>(indexRedI, indexRedJ)
              = BLAIT * hEuk.block<quaternionDim, centerLineDim>(indexI, indexJ);
        }
      }

      for (auto i = 0U; i < feDisp.size(); ++i) {
        const size_t indexRedI = i * centerLineDim;
        const size_t indexI    =  i * centerLineDim;
        for (auto j = 0U; j < feQuats.size(); ++j) {
          const size_t indexRedJ  = dispHessianSize + j * quaternionCorrectionDim;
          const size_t indexJ     = dispHessianSize + j * quaternionDim;
          const auto globalIndexJ = localView_.index(localView_.tree().child(_1, 0).localIndex(j));
          const auto BLAJ         = quaternions[globalIndexJ[1]].orthonormalFrame();

          hred.template block<centerLineDim, quaternionCorrectionDim>(indexRedI, indexRedJ)
              = hEuk.block<centerLineDim, quaternionDim>(indexI, indexJ) * BLAJ;
        }
      }

      hred.topLeftCorner(dispHessianSize,dispHessianSize) = hEuk.topLeftCorner(dispHessianSize,dispHessianSize);
    }

    [[nodiscard]] int size() const { return localViewReduced.size(); }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& rieGrad) const {
      using namespace Dune::Indices;
      Eigen::VectorXdual dx1st;
      dx1st.resize(localView_.size());
      dx1st.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };

      autodiff::dual e;
      autodiff::gradient(f, autodiff::wrt(dx1st), at(dx1st), e, eukGrad);
      const auto& quaternions     = par.getSolution(Ikarus::FESolutions::displacementAndQuaternions)[_1];
      auto& first_child = localView_.tree().child(_0, 0);
      auto& second_child = localView_.tree().child(_1, 0);
      const auto& feDisp = first_child.finiteElement();
      const auto& feQuat = second_child.finiteElement();

      const int dispElementEntries    = feDisp.size() * centerLineDim;
      const int magElementEntries    = feQuat.size() * quaternionDim;
      const int magRedElementEntries = feQuat.size() * quaternionCorrectionDim;
      rieGrad.segment(0,dispElementEntries) = eukGrad.segment(0,dispElementEntries) ;
      for (auto i = 0U; i < feQuat.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        size_t indexRed  = dispElementEntries+ i * quaternionCorrectionDim;
        size_t index     = dispElementEntries +i * quaternionDim;
        rieGrad.template segment<quaternionCorrectionDim>(indexRed)
            = quaternions[globalIndex[1]].orthonormalFrame().transpose() * eukGrad.template segment<quaternionDim>(index);
      }
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(localView_.size());
      dx.setZero();

      return this->calculateScalarImpl(par, dx);
    }

//    void calculateAt(const ResultRequirementsType& req, const Eigen::Vector<double, Traits::mydim>& local,
//                     ResultTypeMap<double>& result) const {
//      using namespace Dune::Indices;
//      const auto& mNodal = req.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_0];
//      const auto& ANodal = req.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_1];
//      const auto& lambda = req.getParameter(Ikarus::FEParameter::loadfactor);
//      auto& child0       = localView_.tree().child(_0, 0);
//      const auto& fe0    = child0.finiteElement();
//      auto& child1       = localView_.tree().child(_1, 0);
//
//      const auto& fe1 = child1.finiteElement();
//      Dune::BlockVector<std::remove_cvref_t<decltype(mNodal[0])>> mN(fe0.size());
//      Dune::BlockVector<std::remove_cvref_t<decltype(ANodal[0])>> AN(fe1.size());
//      for (auto i = 0U; i < fe0.size(); ++i) {
//        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
//        mN[i]            = mNodal[globalIndex[1]];
//      }
//
//      const int magElementEntries = fe0.size() * directorDim;
//      for (auto i = 0U; i < fe1.size(); ++i) {
//        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
//        AN[i]            = ANodal[globalIndex[1]];
//      }
//
//      const auto geo = localView_.element().geometry();
//
//      auto gp = toFieldVector(local);
//
//      using namespace Ikarus::DerivativeDirections;
//      Ikarus::ProjectionBasedLocalFunction magnetLocalFunction(localBasisCenterLine, mN);
//
//      const auto J    = toEigenMatrix(geo.jacobianTransposed(gp)).transpose().eval();
//      const auto Jinv = J.inverse().eval();
//      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisQuaternions, AN);
//
//      const auto normalizedMag = magnetLocalFunction.evaluateFunction(gp).getValue();
//      const auto vectorPot     = vectorPotLocalFunction.evaluateFunction(gp).getValue();
//
//      const Eigen::Matrix<double, directorDim, Traits::mydim> gradm
//          = magnetLocalFunction.evaluateDerivative(gp, wrt(spatialAll), transformWith(Jinv));
//      const Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
//          = vectorPotLocalFunction.evaluateDerivative(gp, wrt(spatialAll), transformWith(Jinv));
//
//      const Eigen::Vector<double, 3> curlA          = Impl::jacobianToCurl(gradA);
//      const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp), lambda);
//
//      typename ResultTypeMap<double>::ResultArray resv;
//      if (req.isResultRequested(ResultType::gradientNormOfMagnetization)) {
//        resv.resize(1, 1);
//        resv(0, 0) = isInside ? gradm.norm() : 0.0;
//        result.insertOrAssignResult(ResultType::gradientNormOfMagnetization, resv);
//      }
//      if (req.isResultRequested(ResultType::BField)) {
//        resv = curlA / (material.mu0 * material.ms) * std::sqrt(2.0);
//        result.insertOrAssignResult(ResultType::BField, resv);
//      }
//      if (req.isResultRequested(ResultType::HField)) {
//        resv.setZero(3, 1);
//        resv(0, 0) = curlA[0] / (Dune::power(material.mu0, 2) * material.ms) * std::sqrt(2.0)
//                     - normalizedMag[0] * material.ms * static_cast<double>(isInside);
//        resv(1, 0) = curlA[1] / (Dune::power(material.mu0, 2) * material.ms) * std::sqrt(2.0)
//                     - normalizedMag[1] * material.ms * static_cast<double>(isInside);
//        if constexpr (directorDim == 3)
//          resv(2, 0) = curlA[2] / (Dune::power(material.mu0, 2) * material.ms) * std::sqrt(2.0)
//                       - normalizedMag[2] * material.ms * static_cast<double>(isInside);
//        result.insertOrAssignResult(ResultType::HField, resv);
//      }
//      if (req.isResultRequested(ResultType::divergenceOfVectorPotential)) {
//        resv.setZero(1, 1);
//        resv(0, 0) = gradA.diagonal().template segment<vectorPotDim>(0).sum();
//        result.insertOrAssignResult(ResultType::divergenceOfVectorPotential, resv);
//      }
//    }

  private:
    template <class ScalarType>
    ScalarType calculateScalarImpl(const FERequirementType& par, Eigen::VectorX<ScalarType>& dx_) const {
      using namespace Dune::Indices;
      const auto& dispNodal = par.getSolution(Ikarus::FESolutions::displacementAndQuaternions)[_0];
      const auto& quatNodal = par.getSolution(Ikarus::FESolutions::displacementAndQuaternions)[_1];
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      auto& child0    = localView_.tree().child(_0, 0);
      const auto& fe0 = child0.finiteElement();
      auto& child1    = localView_.tree().child(_1, 0);
      const auto& fe1 = child1.finiteElement();

      Dune::BlockVector<typename std::remove_cvref_t<decltype(dispNodal[0])>::template Rebind<ScalarType>::other>
          dispEleNodal(
          fe0.size());
      Dune::BlockVector<typename std::remove_cvref_t<decltype(quatNodal[0])>::template Rebind<ScalarType>::other>
          quatEleNodal(
          fe1.size());

      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        dispEleNodal[i].setValue(dispNodal[globalIndex[1]].getValue()+dx_.template segment<centerLineDim>(i * centerLineDim)) ;
      }

      const int dispElementEntries = fe0.size() * centerLineDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        quatEleNodal[i].setValue(quatNodal[globalIndex[1]].getValue()+ dx_.template segment<quaternionDim>(dispElementEntries + i * quaternionDim));
      }

      const auto geo = localView_.element().geometry();

      const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
      Ikarus::StandardLocalFunction centerLineF(localBasisCenterLine, dispEleNodal);
      Ikarus::ProjectionBasedLocalFunction quatF(localBasisQuaternions, quatEleNodal);
      using namespace Ikarus::DerivativeDirections;
      ScalarType energy = 0.0;
      for (const auto& gp : rule) {
        const auto a1            = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        const auto Jinv          = a1.norm();
        const auto r             = centerLineF.evaluateFunction(gp.position());
        const auto rGrad             = centerLineF.evaluateDerivative(gp.position(),wrt(spatialAll), transformWith(Jinv));
        const auto q             = quatF.evaluateFunction(gp.position());
        const auto qE = Eigen::Quaternion<ScalarType>(q);
        const auto qGrad             = quatF.evaluateDerivative(gp.position(),wrt(spatialAll), transformWith(Jinv));

        const auto D = rotationMatrixColumnDerivatives(q);
         Eigen::Matrix<ScalarType,3,3> RPrime;
         RPrime << D[0]*qGrad, D[2]*qGrad, D[2]*qGrad;

         const auto R = qE.toRotationMatrix().eval();
         const auto kappa = (RPrime*R).eval();
         const Eigen::Vector<ScalarType,3> kappaV({kappa(2,1),kappa(0,2),kappa(1,0)});
         const auto eps = (R*rGrad-Eigen::Vector<double,3>::UnitZ()).eval();
          const double G = Ikarus::convertLameConstants({.emodul=material.E,.nu=material.nu}).toShearModulus();
          const Eigen::Vector3d dM({material.E*material.I1,material.E*material.I2,material.J});
          const Eigen::Vector3d cM({G*material.A,G*material.A,material.E*material.A});


        energy += (0.5 * (eps.dot(cM.asDiagonal()*eps)+ kappaV.dot(dM.asDiagonal()*kappaV)))
                  * geo.integrationElement(gp.position()) * gp.weight();  // demag energy
      }
      return energy;
    }

    mutable Eigen::Matrix<double, Eigen::Dynamic, Traits::mydim> dNA;

    LocalViewEmbedded localView_;
    LocalViewReduced localViewReduced;
    Ikarus::LocalBasis<std::remove_cvref_t<
        decltype(std::declval<LocalViewEmbedded>().tree().child(Dune::Indices::_0, 0).finiteElement().localBasis())>>
        localBasisCenterLine;
    Ikarus::LocalBasis<std::remove_cvref_t<
        decltype(std::declval<LocalViewEmbedded>().tree().child(Dune::Indices::_1, 0).finiteElement().localBasis())>>
        localBasisQuaternions;
    std::function<Eigen::Vector<double, centerLineDim>(const Eigen::Vector<double, Traits::mydim>&, const double&)>
        volumeLoad;
    BeamMaterial material;
    unsigned int order{};
  };

}  // namespace Ikarus

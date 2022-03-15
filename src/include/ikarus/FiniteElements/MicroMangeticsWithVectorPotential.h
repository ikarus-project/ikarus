
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
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <concepts>
#include <iostream>
#include <numbers>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include "ikarus/Geometries/SimpleLocalFunction.h"
#include "ikarus/LocalBasis/localBasis.h"
#include "ikarus/utils/LinearAlgebraHelper.h"
#include <ikarus/FiniteElements/AutodiffFE.h>
#include <ikarus/FiniteElements/FEPolicies.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/FiniteElements/physicsHelper.h>
#include <ikarus/Geometries/GeometryWithExternalInput.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::FiniteElements {

  namespace Impl {
    template <typename ScalarType, int FieldSize, int worlddim>
    Eigen::Vector<ScalarType, 3> jacobianToCurl(const Eigen::Matrix<ScalarType, FieldSize, worlddim>& jaco) {
      Eigen::Vector<ScalarType, 3> curl;
      if constexpr (FieldSize == 3 and worlddim == 3) {
        curl[0] = jaco(2, 1) - jaco(1, 2);
        curl[1] = jaco(0, 2) - jaco(2, 0);
        curl[2] = jaco(1, 0) - jaco(0, 1);
      } else if constexpr (FieldSize == 3 and worlddim == 2) {
        curl[0] = jaco(2, 1);
        curl[1] = -jaco(2, 0);
        curl[2] = jaco(1, 0) - jaco(0, 1);
      } else if constexpr (FieldSize == 1 and worlddim == 2) {
        curl[0] = jaco(0, 1);
        curl[1] = -jaco(0, 0);
        curl[2] = 0;
      }
      return curl;
    }
  }  // namespace Impl

  struct MagneticMaterial {
    double A{0.0};
    double K{0.0};
    double mu0{4 * ::std::numbers::pi * 1e-7};
    double ms{0.0};
  };
  template <typename BasisEmbedded, typename BasisReduced>
  class MicroMagneticsWithVectorPotential {
  public:
    static constexpr int directorDim           = BasisEmbedded::PreBasis::template SubPreBasis<0>::Node::CHILDREN;
    static constexpr int vectorPotDim          = BasisEmbedded::PreBasis::template SubPreBasis<1>::Node::CHILDREN;
    static constexpr int directorCorrectionDim = directorDim - 1;
    using DirectorVector                       = Dune::BlockVector<Ikarus::UnitVector<double, directorDim>>;
    using VectorPotVector                      = Dune::BlockVector<Ikarus::RealTuple<double, vectorPotDim>>;
    using MultiTypeVector                      = Dune::MultiTypeBlockVector<DirectorVector, VectorPotVector>;

    using FERequirementType      = FErequirements<MultiTypeVector>;
    using ResultRequirementsType = ResultRequirements<FERequirementType>;
    using LocalViewEmbedded      = typename BasisEmbedded::LocalView;
    using LocalViewReduced       = typename BasisReduced::LocalView;

    template <typename VolumeLoad>
    MicroMagneticsWithVectorPotential(BasisEmbedded& globalBasis, BasisReduced& globalBasisRed,
                                      const typename LocalViewEmbedded::Element& element, MagneticMaterial& p_material,
                                      const VolumeLoad& p_volumeLoad)
        : localView_{globalBasis.localView()},
          localViewReduced{globalBasisRed.localView()},
          volumeLoad(p_volumeLoad),
          material{p_material} {
      localView_.bind(element);
      localViewReduced.bind(element);
      const int order = 3 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();
      //      std::cout<<"Order:"<<order<<std::endl;
      localBasisMag = Ikarus::LocalBasis(localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis());
      localBasisMag.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order), 0, 1);
      localBasisVecPot = Ikarus::LocalBasis(localView_.tree().child(Dune::Indices::_1, 0).finiteElement().localBasis());
      localBasisVecPot.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order), 0,
                            1);
    }

    using Traits = TraitsFromLocalView<LocalViewEmbedded>;
    template <typename ST>
    using LocalFuncMag = Ikarus::Geometry::SimpleLocalFunction<ST, directorDim, Traits::mydim>;
    template <typename ST>
    using LocalFuncVecPot = Ikarus::Geometry::SimpleLocalFunction<ST, vectorPotDim, Traits::mydim>;

    using GlobalIndex = typename LocalViewReduced::MultiIndex;
    void globalIndices(std::vector<GlobalIndex>& globalIndices) const {
      using namespace Dune::Indices;
      const auto& fe = localViewReduced.tree().child(_0, 0).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < directorCorrectionDim; ++j) {
          globalIndices.push_back(localViewReduced.index((localViewReduced.tree().child(_0, j).localIndex(i))));
        }
      }
      const auto& fe2 = localViewReduced.tree().child(_1, 0).finiteElement();
      for (size_t i = 0; i < fe2.size(); ++i) {
        for (int j = 0; j < vectorPotDim; ++j) {
          globalIndices.push_back(localViewReduced.index((localViewReduced.tree().child(_1, j).localIndex(i))));
        }
      }
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& hred) const {
      using namespace Dune::Indices;
      dx2nd.resize(localView_.size());
      dx2nd.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      autodiff::dual2nd e;
      autodiff::hessian(f, wrt(dx2nd), at(dx2nd), e, eukGrad, hEuk);
      Eigen::MatrixXd hessTest;
      hessTest.template resizeLike(hEuk);
      calculateEuclideanHessian(par, hessTest);
      std::cout<<hessTest<<std::endl;
      std::cout<<hEuk<<std::endl;
      calculateEuclideanGradient(par, eukGrad);
      const auto& m        = par.sols[0].get()[_0];
      const auto& feMag    = localView_.tree().child(_0, 0).finiteElement();
      const auto& feVecPot = localView_.tree().child(_1, 0).finiteElement();
      for (auto i = 0U; i < feMag.size(); ++i) {
        const size_t indexRedI  = i * directorCorrectionDim;
        const size_t indexI     = i * directorDim;
        const auto globalIndexI = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        const Eigen::Matrix<double, directorCorrectionDim, directorDim> BLAIT
            = m[globalIndexI[1]].orthonormalFrame().transpose();
        for (auto j = 0U; j < feMag.size(); ++j) {
          const size_t indexRedJ  = j * directorCorrectionDim;
          const size_t indexJ     = j * directorDim;
          const auto globalIndexJ = localView_.index(localView_.tree().child(_0, 0).localIndex(j));
          const auto BLAJ         = m[globalIndexJ[1]].orthonormalFrame();

          hred.template block<directorCorrectionDim, directorCorrectionDim>(indexRedI, indexRedJ)
              = BLAIT * hEuk.block<directorDim, directorDim>(indexI, indexJ) * BLAJ;
        }
        hred.template block<directorCorrectionDim, directorCorrectionDim>(indexRedI, indexRedI)
            -= m[globalIndexI[1]].getValue().dot(eukGrad.template segment<directorDim>(indexI))
               * Eigen::Matrix<double, directorCorrectionDim, directorCorrectionDim>::Identity();
      }

      const int magHessianSize     = feMag.size() * directorCorrectionDim;
      const int magFullHessianSize = feMag.size() * directorDim;

      for (auto i = 0U; i < feMag.size(); ++i) {
        const size_t indexRedI  = i * directorCorrectionDim;
        const size_t indexI     = i * directorDim;
        const auto globalIndexI = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        const Eigen::Matrix<double, directorCorrectionDim, directorDim> BLAIT
            = m[globalIndexI[1]].orthonormalFrame().transpose();
        for (auto j = 0U; j < feVecPot.size(); ++j) {
          const size_t indexRedJ = magHessianSize + j * vectorPotDim;
          const size_t indexJ    = magFullHessianSize + j * vectorPotDim;
          hred.template block<directorCorrectionDim, vectorPotDim>(indexRedI, indexRedJ)
              = BLAIT * hEuk.block<directorDim, vectorPotDim>(indexI, indexJ);
        }
      }

      for (auto i = 0U; i < feVecPot.size(); ++i) {
        const size_t indexRedI = magHessianSize + i * vectorPotDim;
        const size_t indexI    = magFullHessianSize + i * vectorPotDim;
        for (auto j = 0U; j < feMag.size(); ++j) {
          const size_t indexRedJ  = j * directorCorrectionDim;
          const size_t indexJ     = j * directorDim;
          const auto globalIndexJ = localView_.index(localView_.tree().child(_0, 0).localIndex(j));
          const auto BLAJ         = m[globalIndexJ[1]].orthonormalFrame();

          hred.template block<vectorPotDim, directorCorrectionDim>(indexRedI, indexRedJ)
              = hEuk.block<vectorPotDim, directorDim>(indexI, indexJ) * BLAJ;
        }
      }

      hred(Eigen::seq(magHessianSize, Eigen::last), Eigen::seq(magHessianSize, Eigen::last))
          = hEuk(Eigen::seq(magFullHessianSize, Eigen::last), Eigen::seq(magFullHessianSize, Eigen::last));
    }

    [[nodiscard]] int size() const { return localViewReduced.size(); }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& rieGrad) const {
      using namespace Dune::Indices;
      dx1st.resize(localView_.size());
      dx1st.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };

      autodiff::dual e;
      //      autodiff::gradient(f, wrt(dx1st), at(dx1st), e, eukGrad);
      //      std::cout<<eukGrad.transpose()<<std::endl;
      calculateEuclideanGradient(par, eukGrad);
      const auto& m = par.sols[0].get()[_0];
      std::vector<Ikarus::UnitVector<double, directorDim> const*> mLocal;
      auto& first_child = localView_.tree().child(_0, 0);
      const auto& feMag = first_child.finiteElement();

      const int magElementEntries    = feMag.size() * directorDim;
      const int magRedElementEntries = feMag.size() * directorCorrectionDim;
      for (auto i = 0U; i < feMag.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        size_t indexRed  = i * directorCorrectionDim;
        size_t index     = i * directorDim;
        rieGrad.template segment<directorCorrectionDim>(indexRed)
            = m[globalIndex[1]].orthonormalFrame().transpose() * eukGrad.template segment<directorDim>(index);
      }
      rieGrad(Eigen::seq(magRedElementEntries, Eigen::last)) = eukGrad(Eigen::seq(magElementEntries, Eigen::last));
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(localView_.size());
      dx.setZero();

      return this->calculateScalarImpl(par, dx);
    }

    void calculateAt(const ResultRequirementsType& res, const Eigen::Vector<double, Traits::mydim>& local,
                     Eigen::VectorXd& result) const {
      using namespace Dune::Indices;
      const auto& mNodal = res.req.sols[0].get()[_0];
      const auto& ANodal = res.req.sols[0].get()[_1];
      const auto& lambda = res.req.parameter.at(FEParameter::loadfactor);

      auto& child0    = localView_.tree().child(_0, 0);
      const auto& fe0 = child0.finiteElement();
      auto& child1    = localView_.tree().child(_1, 0);

      const auto& fe1 = child1.finiteElement();
      Eigen::Matrix<double, directorDim, Eigen::Dynamic> mN;
      Eigen::Matrix<double, vectorPotDim, Eigen::Dynamic> AN;
      mN.setZero(Eigen::NoChange, fe0.size());
      AN.setZero(Eigen::NoChange, fe1.size());
      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN.col(i)        = mNodal[globalIndex[1]].getValue();
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN.col(i)        = ANodal[globalIndex[1]].getValue();
      }

      const auto geo = localView_.element().geometry();

      const int order = 4 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();

      const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);

      auto gp = toFieldVector(local);
      localBasisMag.template evaluateFunctionAndJacobian(gp, Nm, dNm);
      localBasisVecPot.template evaluateFunctionAndJacobian(gp, NA, dNA);

      const auto J = toEigenMatrix(geo.jacobianTransposed(gp)).transpose().eval();
      Eigen::Vector<double, directorDim> magnetizationEuk;
      Eigen::Vector<double, vectorPotDim> vectorPot;
      magnetizationEuk.setZero();

      for (int i = 0; i < Nm.size(); ++i)
        magnetizationEuk += mN.col(i) * Nm[i];
      const Eigen::Vector<double, directorDim> normalizedMag = magnetizationEuk.normalized();

      vectorPot.setZero();

      for (int i = 0; i < Nm.size(); ++i)
        vectorPot += AN.col(i) * NA[i];

      const double mLength = magnetizationEuk.norm();
      const Eigen::Matrix<double, directorDim, directorDim> Pm
          = (Eigen::Matrix<double, directorDim, directorDim>::Identity() - normalizedMag * normalizedMag.transpose())
            / mLength;
      const auto dNmdx = (dNm * J.inverse()).eval();
      const auto dNAdx = (dNA * J.inverse()).eval();
      Eigen::Matrix<double, directorDim, Traits::mydim> gradm
          = (LocalFuncMag<double>::jacobianTransposed(dNmdx, mN).transpose());
      gradm = Pm * gradm;
      Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
          = (LocalFuncVecPot<double>::jacobianTransposed(dNAdx, AN).transpose());

      const Eigen::Vector<double, 3> curlA(gradA(0, 1), -gradA(0, 0), 0);
      const double divA                             = 0;
      const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp), lambda);
      switch (res.resType) {
        case ResultType::gradientNormOfMagnetization:
          result.resize(1);
          result[0] = (gradm.transpose() * gradm).trace();
          break;
        default:
          DUNE_THROW(Dune::NotImplemented, "This result type is not implemented by this element");
      }
    }

  private:
    void calculateEuclideanGradient(const FERequirementType& par, typename Traits::VectorType& eukGrad_) const {
      eukGrad_.setZero();
      using namespace Dune::Indices;
      const auto& mNodal = par.sols[0].get()[_0];
      const auto& ANodal = par.sols[0].get()[_1];
      const auto& lambda = par.parameter.at(FEParameter::loadfactor);

      auto& child0    = localView_.tree().child(_0, 0);
      const auto& fe0 = child0.finiteElement();
      auto& child1    = localView_.tree().child(_1, 0);

      const auto& fe1 = child1.finiteElement();
      Eigen::Matrix<double, directorDim, Eigen::Dynamic> mN;
      Eigen::Matrix<double, vectorPotDim, Eigen::Dynamic> AN;
      mN.setZero(Eigen::NoChange, fe0.size());
      AN.setZero(Eigen::NoChange, fe1.size());
      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN.col(i)        = mNodal[globalIndex[1]].getValue();
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN.col(i)        = ANodal[globalIndex[1]].getValue();
      }

      const auto geo = localView_.element().geometry();

      const int order = 4 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();

      const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);

      for (const auto& gp : rule) {
        localBasisMag.template evaluateFunctionAndJacobian(gp.position(), Nm, dNm);
        localBasisVecPot.template evaluateFunctionAndJacobian(gp.position(), NA, dNA);

        const auto J = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        Eigen::Vector<double, directorDim> magnetizationEuk;
        Eigen::Vector<double, vectorPotDim> vectorPot;
        magnetizationEuk.setZero();

        for (int i = 0; i < Nm.size(); ++i)
          magnetizationEuk += mN.col(i) * Nm[i];
        const Eigen::Vector<double, directorDim> normalizedMag = magnetizationEuk.normalized();

        vectorPot.setZero();

        for (int i = 0; i < Nm.size(); ++i)
          vectorPot += AN.col(i) * NA[i];

        const double mLength = magnetizationEuk.norm();
        const Eigen::Matrix<double, directorDim, directorDim> Pm
            = (Eigen::Matrix<double, directorDim, directorDim>::Identity() - normalizedMag * normalizedMag.transpose())
              / mLength;
        const auto dNmdx = (dNm * J.inverse()).eval();
        const auto dNAdx = (dNA * J.inverse()).eval();
        Eigen::Matrix<double, directorDim, Traits::mydim> gradw
            = (LocalFuncMag<double>::jacobianTransposed(dNmdx, mN).transpose());
        const Eigen::Matrix<double, directorDim, Traits::mydim> gradm = Pm * gradw;
        Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
            = (LocalFuncVecPot<double>::jacobianTransposed(dNAdx, AN).transpose());

        const Eigen::Vector<double, 3> curlA(gradA(0, 1), -gradA(0, 0), 0);
        const double divA        = 0;
        const double invwsquared = mLength * mLength;
        const Eigen::Matrix<double, directorDim, directorDim> eye
            = Eigen::Matrix<double, directorDim, directorDim>::Identity();
        const Eigen::Matrix<double, directorDim, directorDim> Q1
            = 1 / invwsquared
              * (normalizedMag.dot(gradw.col(0)) * (3 * (normalizedMag * normalizedMag.transpose()) - eye)
                 - gradw.col(0) * normalizedMag.transpose() - normalizedMag * gradw.col(0).transpose());
        const Eigen::Matrix<double, directorDim, directorDim> Q2
            = 1 / invwsquared
              * (normalizedMag.dot(gradw.col(1)) * (3 * (normalizedMag * normalizedMag.transpose()) - eye)
                 - gradw.col(1) * normalizedMag.transpose() - normalizedMag * gradw.col(1).transpose());

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);
        //        std::cout<<"Hbar::"<<Hbar<<std::endl;
        for (size_t i = 0; i < fe0.size(); ++i) {
          const int index                                           = i * directorDim;
          const Eigen::Matrix<double, directorDim, directorDim> WI1 = Q1 * Nm[i] + Pm * dNmdx(i, 0);
          const Eigen::Matrix<double, directorDim, directorDim> WI2 = Q2 * Nm[i] + Pm * dNmdx(i, 1);
          eukGrad_.template segment<directorDim>(index)
              += ((WI1 * gradm.col(0) + WI2 * gradm.col(1) - Pm * curlA * Nm[i]))
                 * geo.integrationElement(gp.position()) * gp.weight();

          eukGrad_.template segment<directorDim>(index)
              -= (Pm * Hbar * Nm[i] / material.ms) * geo.integrationElement(gp.position()) * gp.weight();
        }
        const int magEukSize = fe0.size() * directorDim;
        for (size_t i = 0; i < fe1.size(); ++i) {
          const int index = magEukSize + i * vectorPotDim;
          Eigen::Vector<double, directorDim> gradCurlA_dI;
          if constexpr (directorDim == 3) {
            gradCurlA_dI[0] = dNAdx(i, 1);
            gradCurlA_dI[1] = -dNAdx(i, 0);
            gradCurlA_dI[2] = 0;
          } else if constexpr (directorDim == 2) {
            gradCurlA_dI[0] = dNAdx(i, 1);
            gradCurlA_dI[1] = -dNAdx(i, 0);
          }
          eukGrad_.template segment<vectorPotDim>(index).array()
              += (-gradCurlA_dI.dot(normalizedMag) + gradCurlA_dI.dot(curlA)) * geo.integrationElement(gp.position())
                 * gp.weight();
        }
      }
      std::cout << eukGrad_.transpose() << std::endl;
    }

    void calculateEuclideanHessian(const FERequirementType& par, typename Traits::MatrixType& eukHess_) const {
      eukHess_.setZero();
      using namespace Dune::Indices;
      const auto& mNodal = par.sols[0].get()[_0];
      const auto& ANodal = par.sols[0].get()[_1];
      const auto& lambda = par.parameter.at(FEParameter::loadfactor);

      auto& child0    = localView_.tree().child(_0, 0);
      const auto& fe0 = child0.finiteElement();
      auto& child1    = localView_.tree().child(_1, 0);

      const auto& fe1 = child1.finiteElement();
      Eigen::Matrix<double, directorDim, Eigen::Dynamic> mN;
      Eigen::Matrix<double, vectorPotDim, Eigen::Dynamic> AN;
      mN.setZero(Eigen::NoChange, fe0.size());
      AN.setZero(Eigen::NoChange, fe1.size());
      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN.col(i)        = mNodal[globalIndex[1]].getValue();
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN.col(i)        = ANodal[globalIndex[1]].getValue();
      }

      const auto geo = localView_.element().geometry();

      const int order = 4 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();

      const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);

      for (const auto& gp : rule) {
        localBasisMag.template evaluateFunctionAndJacobian(gp.position(), Nm, dNm);
        localBasisVecPot.template evaluateFunctionAndJacobian(gp.position(), NA, dNA);

        const auto J = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        Eigen::Vector<double, directorDim> magnetizationEuk;
        Eigen::Vector<double, vectorPotDim> vectorPot;
        magnetizationEuk.setZero();

        for (int i = 0; i < Nm.size(); ++i)
          magnetizationEuk += mN.col(i) * Nm[i];
        const Eigen::Vector<double, directorDim> normalizedMag = magnetizationEuk.normalized();

        vectorPot.setZero();

        for (int i = 0; i < Nm.size(); ++i)
          vectorPot += AN.col(i) * NA[i];

        const double mLength = magnetizationEuk.norm();
        const Eigen::Matrix<double, directorDim, directorDim> Pm
            = (Eigen::Matrix<double, directorDim, directorDim>::Identity() - normalizedMag * normalizedMag.transpose())
              / mLength;
        const auto dNmdx = (dNm * J.inverse()).eval();
        const auto dNAdx = (dNA * J.inverse()).eval();
        Eigen::Matrix<double, directorDim, Traits::mydim> gradw
            = (LocalFuncMag<double>::jacobianTransposed(dNmdx, mN).transpose());
        const Eigen::Matrix<double, directorDim, Traits::mydim> gradm = Pm * gradw;
        Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
            = (LocalFuncVecPot<double>::jacobianTransposed(dNAdx, AN).transpose());

        const Eigen::Vector<double, 3> curlA(gradA(0, 1), -gradA(0, 0), 0);
        const double divA           = 0;
        const double invLength = 1/mLength;
        const double invwsquared    = invLength * invLength;
        using DirectorMatrix        = Eigen::Matrix<double, directorDim, directorDim>;
        const DirectorMatrix eye    = DirectorMatrix::Identity();
        const DirectorMatrix mdyadm = normalizedMag * normalizedMag.transpose();
        const DirectorMatrix Q1
            =  invwsquared
              * (normalizedMag.dot(gradw.col(0)) * (3 * (mdyadm)-eye) - gradw.col(0) * normalizedMag.transpose()
                 - normalizedMag * gradw.col(0).transpose());
        const DirectorMatrix Q2
            =  invwsquared
              * (normalizedMag.dot(gradw.col(1)) * (3 * (mdyadm)-eye) - gradw.col(1) * normalizedMag.transpose()
                 - normalizedMag * gradw.col(1).transpose());

        const DirectorMatrix Id3minus5tdyadt = eye - 5.0 * mdyadm;
        const double normwcubinv             =  invwsquared * invLength;
        const DirectorMatrix td1dyadt        = gradm.col(0) * normalizedMag.transpose();
        const DirectorMatrix td2dyadt        = gradm.col(1) * normalizedMag.transpose();

        const double td1scalwd1 = gradm.col(0).dot(gradw.col(0));
        const double td1scalwd2 = gradm.col(0).dot(gradw.col(1));
        const double td2scalwd1 = gradm.col(1).dot(gradw.col(0));
        const double td2scalwd2 = gradm.col(1).dot(gradw.col(1));
        const double tscalwd1   = normalizedMag.dot(gradw.col(0));
        const double tscalwd2   = normalizedMag.dot(gradw.col(1));
        DirectorMatrix chi11d   = normwcubinv
                                * (3 * tscalwd1 * td1dyadt + 3 * (0.5 * td1scalwd1 * mdyadm)
                                   - gradm.col(0) * gradw.col(0).transpose() - td1scalwd1 * 0.5 * eye);
        chi11d                = (chi11d + chi11d.transpose()).eval();
        DirectorMatrix chi22d = normwcubinv
                                * (3 * tscalwd2 * td2dyadt + 3 * (0.5 * td2scalwd2 * mdyadm)
                                   - gradm.col(1) * gradw.col(1).transpose() - td2scalwd2 * 0.5 * eye);
        chi22d                = (chi22d + chi22d.transpose()).eval();
        DirectorMatrix chi12d = normwcubinv
                                * (3 * tscalwd2 * td1dyadt + 3 * (0.5 * td1scalwd2 * mdyadm)
                                   - gradm.col(0) * gradw.col(1).transpose() - td1scalwd2 * 0.5 * eye);
        chi12d                = (chi12d + chi12d.transpose()).eval();
        DirectorMatrix chi21d = normwcubinv
                                * (3 * tscalwd1 * td2dyadt + 3 * (0.5 * td2scalwd1 * mdyadm)
                                   - gradm.col(1) * gradw.col(0).transpose() - td2scalwd1 * 0.5 * eye);
        chi21d = (chi21d + chi21d.transpose()).eval();

        const DirectorMatrix S1d    =  invwsquared * (-td1dyadt - (normalizedMag * gradm.col(0).transpose()));
        const DirectorMatrix S2d    =  invwsquared * (-td2dyadt - (normalizedMag * gradm.col(1).transpose()));
        const DirectorMatrix chifac = chi11d + chi22d + chi12d + chi21d;

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);
        //        std::cout<<"Hbar::"<<Hbar<<std::endl;
        for (size_t i = 0; i < fe0.size(); ++i) {
          const int indexI                                          = i * directorDim;
          const Eigen::Matrix<double, directorDim, directorDim> WI1 = Q1 * Nm[i] + Pm * dNmdx(i, 0);
          const Eigen::Matrix<double, directorDim, directorDim> WI2 = Q2 * Nm[i] + Pm * dNmdx(i, 1);
          for (size_t j = 0; j < fe0.size(); ++j) {
            const int indexJ                                          = j * directorDim;
            const Eigen::Matrix<double, directorDim, directorDim> WJ1 = Q1 * Nm[j] + Pm * dNmdx(j, 0);
            const Eigen::Matrix<double, directorDim, directorDim> WJ2 = Q2 * Nm[j] + Pm * dNmdx(j, 1);
            const double NdN1                                         = dNmdx(j, 0) * Nm[i] + Nm[j] * dNmdx(i, 0);
            const double NdN2                                         = dNmdx(j, 1) * Nm[i] + Nm[j] * dNmdx(i, 1);
            eukHess_.template block<directorDim, directorDim>(indexI, indexJ)
                += WI1 * WJ1 + WI2 * WJ2 + (WI1 * WJ2 + WI2 * WJ1) + Nm[i] * Nm[j] * chifac
                   + (S1d * NdN1 + S2d * NdN2 + (S1d * NdN2 + S2d * NdN1));

          }
        }
//        const int magEukSize = fe0.size() * directorDim;
//        for (size_t i = 0; i < fe1.size(); ++i) {
//          const int index = magEukSize + i * vectorPotDim;
//          Eigen::Vector<double, directorDim> gradCurlA_dI;
//          if constexpr (directorDim == 3) {
//            gradCurlA_dI[0] = dNAdx(i, 1);
//            gradCurlA_dI[1] = -dNAdx(i, 0);
//            gradCurlA_dI[2] = 0;
//          } else if constexpr (directorDim == 2) {
//            gradCurlA_dI[0] = dNAdx(i, 1);
//            gradCurlA_dI[1] = -dNAdx(i, 0);
//          }
//          eukGrad_.template segment<vectorPotDim>(index).array()
//              += (-gradCurlA_dI.dot(normalizedMag) + gradCurlA_dI.dot(curlA)) * geo.integrationElement(gp.position())
//                 * gp.weight();
//        }
      }

    }

    template <class ScalarType>
    ScalarType calculateScalarImpl(const FERequirementType& par, Eigen::VectorX<ScalarType>& dx_) const {
      using namespace Dune::Indices;
      const auto& mNodal = par.sols[0].get()[_0];
      const auto& ANodal = par.sols[0].get()[_1];
      const auto& lambda = par.parameter.at(FEParameter::loadfactor);

      auto& child0    = localView_.tree().child(_0, 0);
      const auto& fe0 = child0.finiteElement();
      auto& child1    = localView_.tree().child(_1, 0);

      const auto& fe1 = child1.finiteElement();
      Eigen::Matrix<ScalarType, directorDim, Eigen::Dynamic> mN;
      Eigen::Matrix<ScalarType, vectorPotDim, Eigen::Dynamic> AN;
      mN.setZero(Eigen::NoChange, fe0.size());
      AN.setZero(Eigen::NoChange, fe1.size());
      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN.col(i)        = dx_.template segment<directorDim>(i * directorDim) + mNodal[globalIndex[1]].getValue();
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN.col(i)        = dx_.template segment<vectorPotDim>(magElementEntries + i * vectorPotDim)
                    + ANodal[globalIndex[1]].getValue();
      }

      ScalarType energy = 0.0;

      const auto geo = localView_.element().geometry();

      const int order = 4 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();

      const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);

      for (const auto& gp : rule) {
        localBasisMag.template evaluateFunctionAndJacobian(gp.position(), Nm, dNm);
        localBasisVecPot.template evaluateFunctionAndJacobian(gp.position(), NA, dNA);

        const auto J = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        Eigen::Vector<ScalarType, directorDim> magnetizationEuk;
        Eigen::Vector<ScalarType, vectorPotDim> vectorPot;
        magnetizationEuk.setZero();

        for (int i = 0; i < Nm.size(); ++i)
          magnetizationEuk += mN.col(i) * Nm[i];
        const Eigen::Vector<ScalarType, directorDim> normalizedMag = magnetizationEuk.normalized();

        vectorPot.setZero();

        for (int i = 0; i < Nm.size(); ++i)
          vectorPot += AN.col(i) * NA[i];

        const ScalarType mLength = magnetizationEuk.norm();
        const Eigen::Matrix<ScalarType, directorDim, directorDim> Pm
            = (Eigen::Matrix<ScalarType, directorDim, directorDim>::Identity()
               - normalizedMag * normalizedMag.transpose())
              / mLength;
        const auto dNmdx = (dNm * J.inverse()).eval();
        const auto dNAdx = (dNA * J.inverse()).eval();
        Eigen::Matrix<ScalarType, directorDim, Traits::mydim> gradm
            = (LocalFuncMag<ScalarType>::jacobianTransposed(dNmdx, mN).transpose());
        gradm = Pm * gradm;
        Eigen::Matrix<ScalarType, vectorPotDim, Traits::mydim> gradA
            = (LocalFuncVecPot<ScalarType>::jacobianTransposed(dNAdx, AN).transpose());

        const Eigen::Vector<ScalarType, 3> curlA(gradA(0, 1), -gradA(0, 0), 0);
        const ScalarType divA = 0;
        //        std::cout<<"mN::"<<mN<<std::endl;
        //        std::cout<<"AN::"<<AN<<std::endl;
        //        std::cout<<"normalizedMag::"<<normalizedMag.transpose()<<std::endl;
        //        std::cout<<"magnetizationEuk::"<<magnetizationEuk.transpose()<<std::endl;
        //        std::cout<<"Pm::"<<Pm<<std::endl;
        //        std::cout<<"gradm::"<<gradm.transpose()<<std::endl;
        //        std::cout<<"dNm::"<<dNm<<std::endl;
        //        std::cout<<"dNmdx::"<<dNmdx<<std::endl;
        //        std::cout<<"gradA::"<<gradA.transpose()<<std::endl;
        //                std::cout<<"curlA::"<<curlA.transpose()<<std::endl;
        //        std::cout<<"divA::"<<divA<<std::endl;
        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);
        //        std::cout<<"Hbar::"<<Hbar<<std::endl;
        energy += (0.5 * (gradm.transpose() * gradm).trace() - 0*2 * normalizedMag.dot(Hbar) / material.ms)
                  * geo.integrationElement(gp.position()) * gp.weight();  // exchange and zeeman energy

        energy += (0.5 * curlA.squaredNorm() - 0*normalizedMag.dot(curlA) + divA * divA)
                  * geo.integrationElement(gp.position()) * gp.weight();  // demag energy
      }
      return energy;
    }

    mutable Eigen::MatrixXd hEuk;
    mutable Eigen::VectorXd eukGrad;
    mutable Eigen::VectorXdual2nd dx2nd;
    mutable Eigen::VectorXdual dx1st;
    mutable Eigen::VectorXd Nm;
    mutable Eigen::Matrix<double, Eigen::Dynamic, Traits::mydim> dNm;
    mutable Eigen::VectorXd NA;
    mutable Eigen::Matrix<double, Eigen::Dynamic, Traits::mydim> dNA;

    LocalViewEmbedded localView_;
    LocalViewReduced localViewReduced;
    Ikarus::LocalBasis<std::remove_cvref_t<
        decltype(std::declval<LocalViewEmbedded>().tree().child(Dune::Indices::_0, 0).finiteElement().localBasis())>>
        localBasisMag;
    Ikarus::LocalBasis<std::remove_cvref_t<
        decltype(std::declval<LocalViewEmbedded>().tree().child(Dune::Indices::_1, 0).finiteElement().localBasis())>>
        localBasisVecPot;
    std::function<Eigen::Vector<double, directorDim>(const Eigen::Vector<double, Traits::mydim>&, const double&)>
        volumeLoad;
    MagneticMaterial material;
  };

}  // namespace Ikarus::FiniteElements

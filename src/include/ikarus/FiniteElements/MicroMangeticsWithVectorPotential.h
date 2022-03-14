
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

    using FERequirementType = FErequirements<MultiTypeVector>;
    using LocalViewEmbedded = typename BasisEmbedded::LocalView;
    using LocalViewReduced  = typename BasisReduced::LocalView;

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
      hessian(f, wrt(dx2nd), at(dx2nd), e, eukGrad, hEuk);

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

      const int magHessianSize = feMag.size() * directorCorrectionDim;
      const int magFullHessianSize = feMag.size() * directorDim;

      for (auto i = 0U; i < feMag.size(); ++i) {
        const size_t indexRedI  = i * directorCorrectionDim;
        const size_t indexI     = i * directorDim;
        const auto globalIndexI = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        const Eigen::Matrix<double, directorCorrectionDim, directorDim> BLAIT
            = m[globalIndexI[1]].orthonormalFrame().transpose();
        for (auto j = 0U; j < feVecPot.size(); ++j) {
          const size_t indexRedJ = magHessianSize + j * vectorPotDim;
          const size_t indexJ = magFullHessianSize + j * vectorPotDim;
          hred.template block<directorCorrectionDim, vectorPotDim>(indexRedI, indexRedJ)
              = BLAIT * hEuk.block<directorDim, vectorPotDim>(indexI, indexJ);
        }
      }

      for (auto i = 0U; i < feVecPot.size(); ++i) {
        const size_t indexRedI = magHessianSize + i * vectorPotDim;
        const size_t indexI = magFullHessianSize + i * vectorPotDim;
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
      autodiff::gradient(f, wrt(dx1st), at(dx1st), e, eukGrad);
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

  private:
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

      const int order = 3 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();

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

        const Eigen::Vector<ScalarType, 3> curlA(gradA(0,1),-gradA(0,0),0);
        const ScalarType divA                    = 0;
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
        energy += (0.5 * (gradm.transpose() * gradm).trace() - 2 * normalizedMag.dot(Hbar) / material.ms)
                  * geo.integrationElement(gp.position()) * gp.weight();  // exchange and zeeman energy

        energy += (0.5 * curlA.squaredNorm() - normalizedMag.dot(curlA) + divA * divA)
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

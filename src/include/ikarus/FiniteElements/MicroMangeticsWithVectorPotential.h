
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

  struct MagneticMaterial {
    double K{0.0};
    double mu0{4 * ::std::numbers::pi * 1e-7};
    double ms{0.0};
  };
  template <typename BasisEmbedded, typename BasisReduced>
  class MicroMagneticsWithVectorPotential {
  public:
    static constexpr int directorDim           = BasisEmbedded::PreBasis::Node::CHILDREN;
    static constexpr int directorCorrectionDim = directorDim - 1;
    using DirectorVector                       = Dune::BlockVector<Ikarus::UnitVector<double, directorDim>>;

    using FERequirementType = FErequirements<DirectorVector>;
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
      const int order = 2;
      localBasis      = Ikarus::LocalBasis(localView_.tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order), 0, 1);
    }

    using Traits = TraitsFromLocalView<LocalViewEmbedded>;
    template <typename ST>
    using DefoGeo = Ikarus::Geometry::GeometryWithExternalInput<ST, directorDim, Traits::mydim>;

    using GlobalIndex = typename LocalViewReduced::MultiIndex;
    [[nodiscard]] std::vector<GlobalIndex> globalIndices() const {
      const auto& fe = localViewReduced.tree().child(0).finiteElement();
      std::vector<GlobalIndex> globalIndices;
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < directorCorrectionDim; ++j) {
          globalIndices.push_back(localViewReduced.index((localViewReduced.tree().child(j).localIndex(i))));
        }
      }
      return globalIndices;
    }

    [[nodiscard]] typename Traits::MatrixType calculateMatrix(const FERequirementType& par) const {
      Eigen::VectorXdual2nd dx(localView_.size());
      dx.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      Eigen::MatrixXd h;
      Eigen::VectorXd g;
      autodiff::dual2nd e;
      hessian(f, wrt(dx), at(dx), e, g, h);
      Eigen::MatrixXd rieHess(localViewReduced.size(), localViewReduced.size());
      const auto& m = par.sols[0].get();
      std::vector<Ikarus::UnitVector<double, directorDim> const*> mLocal;
      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      for (auto i = 0U; i < fe.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(0).localIndex(i));
        mLocal.push_back( &m[globalIndex[0]]);
//        std::cout<<*mLocal.back()<<std::endl;
      }
      for (int i = 0; auto& mSingleI : mLocal) {
        size_t indexRedI = i * directorCorrectionDim;
        size_t indexI    = i * directorDim;
        auto BLAIT       = mSingleI->orthonormalFrame().transpose();
        for (int j = 0; auto& mSingleJ : mLocal) {
          size_t indexRedJ = j * directorCorrectionDim;
          size_t indexJ    = j * directorDim;
          auto BLAJ        = mSingleJ->orthonormalFrame();
          rieHess.block<directorCorrectionDim, directorCorrectionDim>(indexRedI, indexRedJ)
              = BLAIT * h.block<directorDim, directorDim>(indexI, indexJ) * BLAJ;
          ++j;
        }
        rieHess.block<directorCorrectionDim, directorCorrectionDim>(indexRedI, indexRedI)
            -= mSingleI->getValue().dot(g.template segment<directorDim>(indexI))
               * Eigen::Matrix<double, directorCorrectionDim, directorCorrectionDim>::Identity();
        ++i;
      }
      return rieHess;
    }

    [[nodiscard]] typename Traits::VectorType calculateVector(const FERequirementType& par) const {
      Eigen::VectorXdual dx(localView_.size());
      dx.setZero();
      auto f       = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      Eigen::VectorXd eukGrad;
      autodiff::dual e;
      autodiff::gradient(f, wrt(dx), at(dx),e,eukGrad);
      Eigen::VectorXd rieGrad(localViewReduced.size());
      const auto& m = par.sols[0].get();
      std::vector<Ikarus::UnitVector<double, directorDim> const*> mLocal;
      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      for (auto i = 0U; i < fe.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(0).localIndex(i));
        mLocal.push_back( &m[globalIndex[0]]);
      }

      for (int i = 0; auto& mSingle : mLocal) {
        size_t indexRed = i * directorCorrectionDim;
        size_t index    = i * directorDim;
        rieGrad.template segment<directorCorrectionDim>(indexRed)
            = mSingle->orthonormalFrame().transpose() * eukGrad.template segment<directorDim>(index);
      }
      return rieGrad;
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(localView_.size());
      dx.setZero();

      return this->calculateScalarImpl(par, dx);
    }

  private:
    template <class ScalarType>
    ScalarType calculateScalarImpl(const FERequirementType& par, Eigen::VectorX<ScalarType>& dx) const {
      const auto& m      = par.sols[0].get();
      const auto& lambda = par.parameter.at(FEParameter::loadfactor);

      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Eigen::Matrix<ScalarType, directorDim, Eigen::Dynamic> mN;
      mN.setZero(Eigen::NoChange, fe.size());
      for (auto i = 0U; i < fe.size(); ++i)
      {  auto globalIndex = localView_.index(localView_.tree().child(0).localIndex(i));
          //          std::cout<<globalIndex[0]<<" "<<globalIndex[1]<<" "<<"i: "<< i <<" k2: "<<k2<<std::endl;
          mN.col(i) = dx.template segment<directorDim>(i*directorDim) + m[globalIndex[0]].getValue();
        }

//      std::cout << "mN" << std::endl;
//      std::cout << mN << std::endl;
      ScalarType energy = 0.0;
      const int order   = 2 * (fe.localBasis().order());
//      const auto& rule  = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);

      const auto geo = localView_.element().geometry();

      for (const auto& [gpIndex, gp, N, dN] : localBasis.viewOverFunctionAndJacobian()) {
        const auto J = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        Eigen::Vector<ScalarType, directorDim> mV;
        mV.setZero();

        for (int i = 0; i < N.size(); ++i)
          mV += mN.col(i) * N[i];
         Eigen::Vector<ScalarType, directorDim> mVn = mV;
         mVn.normalize();
//        std::cout<<mVn<<std::endl;
        const ScalarType mLength                         = mV.norm();
        const Eigen::Matrix<ScalarType, directorDim, directorDim> Pm
            = (Eigen::Matrix<ScalarType, directorDim, directorDim>::Identity() - mVn * mVn.transpose()) / mLength;
        const auto dNdx  = (dN * J.inverse()).eval();
        const auto gradm = (Pm * DefoGeo<ScalarType>::jacobianTransposed(dNdx, mN).transpose()).eval();

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);
        energy += (0.5 * (gradm.transpose() * gradm).trace() - 0*mVn.dot(Hbar)) * geo.integrationElement(gp.position())
                  * gp.weight();
      }
      return energy;
    }

    LocalViewEmbedded localView_;
    LocalViewReduced localViewReduced;
    Ikarus::LocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalViewEmbedded>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::function<Eigen::Vector<double, directorDim>(const Eigen::Vector<double, Traits::mydim>&,
                                                          const double&)>
        volumeLoad;
    MagneticMaterial material;
  };

}  // namespace Ikarus::FiniteElements

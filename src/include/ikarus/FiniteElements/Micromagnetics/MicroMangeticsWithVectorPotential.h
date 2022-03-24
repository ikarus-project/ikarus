
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

#include "ikarus/FiniteElements/Interface/FEPolicies.h"
#include "ikarus/FiniteElements/Interface/FiniteElementFunctionConcepts.h"
#include "ikarus/FiniteElements/Interface/InterfaceFiniteElement.h"
#include "ikarus/FiniteElements/physicsHelper.h"
#include "ikarus/LocalBasis/localBasis.h"

#include "ikarus/LocalFunctions/ProjectionBasedLocalFunction.h"
#include "ikarus/LocalFunctions/StandardLocalFunction.h"
#include "ikarus/Variables/VariableDefinitions.h"
#include "ikarus/utils/LinearAlgebraHelper.h"
#include "ikarus/utils/LinearAlgebraTypedefs.h"

namespace Ikarus::FiniteElements {

  namespace Impl {
    template <typename Derived>
    Eigen::Vector<typename Derived::Scalar, 3> jacobianToCurl(const Eigen::MatrixBase<Derived>& jaco) {
      Eigen::Vector<typename Derived::Scalar, 3> curl;
      constexpr int FieldSize = Derived::RowsAtCompileTime;
      constexpr int worlddim = Derived::ColsAtCompileTime;
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
      order         = 4 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();
      localBasisMag = Ikarus::LocalBasis(localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis());
      localBasisMag.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order), bindDerivatives(0, 1));
      localBasisVecPot = Ikarus::LocalBasis(localView_.tree().child(Dune::Indices::_1, 0).finiteElement().localBasis());
      localBasisVecPot.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order), bindDerivatives(0, 1));
    }

    using Traits = TraitsFromLocalView<LocalViewEmbedded>;

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
      hEuk.resize(localView_.size(), localView_.size());
      eukGrad.resize(localView_.size());
      calculateEuclideanHessian(par, hEuk);
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
      calculateEuclideanGradient(par, eukGrad);
      const auto& m = par.sols[0].get()[_0];
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
      Dune::BlockVector<std::remove_cvref_t<decltype(mNodal[0])>> mN(fe0.size());
      Dune::BlockVector<std::remove_cvref_t<decltype(ANodal[0])>> AN(fe1.size());
      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN[i]            = mNodal[globalIndex[1]];
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN[i]        = ANodal[globalIndex[1]];
      }

      const auto geo = localView_.element().geometry();

      auto gp = toFieldVector(local);
      localBasisMag.template evaluateFunctionAndJacobian(gp, Nm, dNm);
      localBasisVecPot.template evaluateFunctionAndJacobian(gp, NA, dNA);

      using namespace Ikarus::DerivativeDirections;
      Ikarus::ProjectionBasedLocalFunction magnetLocalFunction(localBasisMag, mN);

      const auto J = toEigenMatrix(geo.jacobianTransposed(gp)).transpose().eval();
      const auto Jinv = J.inverse().eval();
      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisVecPot, AN);

      const auto normalizedMag = magnetLocalFunction.evaluateFunction(gp).getValue();
      const auto vectorPot = vectorPotLocalFunction.evaluateFunction(gp).getValue();

      const Eigen::Matrix<double, directorDim, Traits::mydim> gradm
          = magnetLocalFunction.evaluateDerivative(gp, wrt(spatialall), transformWith(Jinv));
      const Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
          = vectorPotLocalFunction.evaluateDerivative(gp, wrt(spatialall), transformWith(Jinv));

      const Eigen::Vector<double, 3> curlA(gradA(0, 1), -gradA(0, 0), 0);
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
      Dune::BlockVector<std::remove_cvref_t<decltype(mNodal[0])>> mN(fe0.size());
      Dune::BlockVector<std::remove_cvref_t<decltype(ANodal[0])>> AN(fe1.size());
      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN[i]            = mNodal[globalIndex[1]];
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN[i]        = ANodal[globalIndex[1]];
      }

      const auto geo = localView_.element().geometry();
      using namespace Ikarus::DerivativeDirections;
      Ikarus::ProjectionBasedLocalFunction magnetLocalFunction(localBasisMag, mN);
      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisVecPot, AN);

      for (const auto& [gpIndex, gp] : magnetLocalFunction.viewOverIntegrationPoints()) {
        const double weight = geo.integrationElement(gp.position()) * gp.weight();
        localBasisVecPot.evaluateFunctionAndJacobian(gp.position(), NA, dNA);
        const auto J             = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        const auto Jinv             = J.inverse().eval();
        const auto normalizedMag = magnetLocalFunction.evaluateFunction(gpIndex).getValue();
        const auto vectorPot = vectorPotLocalFunction.evaluateFunction(gpIndex).getValue();

        const auto dNAdx = (dNA * Jinv).eval();
        const Eigen::Matrix<double, directorDim, Traits::mydim> gradm
            = magnetLocalFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));
        const Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
            = vectorPotLocalFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));

        const Eigen::Vector<double, 3> curlA=Impl::jacobianToCurl(gradA);

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);

        for (size_t i = 0; i < fe0.size(); ++i) {
          const int index = i * directorDim;
          const auto WI   = magnetLocalFunction.evaluateDerivative(gpIndex, wrt(spatialall, coeffs),
                                                                   transformWith(Jinv), coeffIndices(i));
          const auto PmI  = magnetLocalFunction.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
          Eigen::Vector<double, directorDim> tmp;
          tmp.setZero();
          for (int j = 0; j < Traits ::mydim; ++j)
            tmp += WI[j] * gradm.col(j);

          eukGrad_.template segment<directorDim>(index) += (tmp - PmI * curlA) * weight;
          eukGrad_.template segment<directorDim>(index) -= (2 * PmI * Hbar / material.ms) * weight;
        }
        const int magEukSize = fe0.size() * directorDim;
        for (size_t i = 0; i < fe1.size(); ++i) {
          const int index = magEukSize + i * vectorPotDim;
          Eigen::Vector<double, directorDim> gradCurlA_dI;
          if constexpr (directorDim == 3) {
            gradCurlA_dI<<dNAdx(i, 1), -dNAdx(i, 0), 0.0;
          } else if constexpr (directorDim == 2) {
            gradCurlA_dI<<dNAdx(i, 1), -dNAdx(i, 0);
          }
          eukGrad_.template segment<vectorPotDim>(index).array()
              += (-gradCurlA_dI.dot(normalizedMag) + gradCurlA_dI.dot(curlA)) * weight;
        }
      }
    }

    void calculateEuclideanHessian(const FERequirementType& par, typename Traits::MatrixType& eukHess_) const {
      eukHess_.setZero();

      using namespace Dune::Indices;
      const auto& mNodal = par.sols[0].get()[_0];
      const auto& ANodal = par.sols[0].get()[_1];
      const auto& lambda = par.parameter.at(FEParameter::loadfactor);

      const auto& child0    = localView_.tree().child(_0, 0);
      const auto& fe0 = child0.finiteElement();
      const auto& child1    = localView_.tree().child(_1, 0);
      const auto& fe1 = child1.finiteElement();

      Dune::BlockVector<std::remove_cvref_t<decltype(mNodal[0])>> mN(fe0.size());
      Dune::BlockVector<std::remove_cvref_t<decltype(ANodal[0])>> AN(fe1.size());
      for (auto i = 0U; i < fe0.size(); ++i) {
        const auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN[i]            = mNodal[globalIndex[1]];
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        const auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN[i]       = ANodal[globalIndex[1]];
      }

      const auto geo = localView_.element().geometry();

      using namespace Ikarus::DerivativeDirections;
      Ikarus::ProjectionBasedLocalFunction mLocalF(localBasisMag, mN);
      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisVecPot, AN);

      for (const auto& [gpIndex, gp] : mLocalF.viewOverIntegrationPoints()) {
        const double weight = geo.integrationElement(gp.position()) * gp.weight();
        localBasisVecPot.template evaluateFunctionAndJacobian(gp.position(), NA, dNA);

        const auto J             = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        const auto Jinv          = J.inverse().eval();
        const auto normalizedMag = mLocalF.evaluateFunction(gpIndex).getValue();
        const auto gradm         = mLocalF.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));
        const auto vectorPot =vectorPotLocalFunction.evaluateFunction(gpIndex).getValue();

        const auto dNAdx = (dNA * J.inverse()).eval();
        Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
            = vectorPotLocalFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));

        const Eigen::Vector<double, 3> curlA= Impl::jacobianToCurl(gradA);

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);
        for (size_t i = 0; i < fe0.size(); ++i) {
          const int indexI = i * directorDim;
          const auto WI
              = mLocalF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(i));
          for (size_t j = 0; j < fe0.size(); ++j) {
            const int indexJ = j * directorDim;
            const auto WJ
                = mLocalF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(j));

            const auto mddHbar
                = mLocalF.evaluateDerivative(gpIndex, wrt(coeffs, coeffs), along(2*Hbar/ material.ms+curlA), coeffIndices(i, j));

            Eigen::Matrix<double,directorDim,directorDim> tmp;
            tmp.setZero();
            for (int k = 0; k < Traits::mydim; ++k)
              tmp+=WI[k]*WJ[k] +mLocalF.evaluateDerivative(gpIndex, wrt(spatial(k), coeffs, coeffs), along(gradm.col(k)),
                                                                transformWith(Jinv), coeffIndices(i, j));;

            eukHess_.template block<directorDim, directorDim>(indexI, indexJ) += tmp * weight;
            eukHess_.template block<directorDim, directorDim>(indexI, indexJ) -=  mddHbar* weight;
          }
        }
        const int magEukSize = fe0.size() * directorDim;
        for (size_t i = 0; i < fe1.size(); ++i) {
          const int indexI = magEukSize + i * vectorPotDim;
          const Eigen::Vector<double, directorDim> gradCurlA_dI= Impl::jacobianToCurl(dNAdx.row(i));
          for (size_t j = 0; j < fe1.size(); ++j) {
            const int indexJ = magEukSize + j * vectorPotDim;
            const Eigen::Vector<double, directorDim> gradCurlA_dJ= Impl::jacobianToCurl(dNAdx.row(j));
            eukHess_.template block<vectorPotDim, vectorPotDim>(indexI, indexJ)
                += (gradCurlA_dI.transpose() * gradCurlA_dJ) * geo.integrationElement(gp.position()) * gp.weight();
          }
        }

        for (size_t i = 0; i < fe0.size(); ++i) {
          const int indexI = i * directorDim;
          const auto PmI   = mLocalF.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
          for (size_t j = 0; j < fe1.size(); ++j) {
            const int indexJ = magEukSize + j * vectorPotDim;
            const Eigen::Vector<double, directorDim> gradCurlA_dJ= Impl::jacobianToCurl(dNAdx.row(j));

            eukHess_.template block<directorDim, vectorPotDim>(indexI, indexJ)
                -= (PmI * gradCurlA_dJ) * geo.integrationElement(gp.position()) * gp.weight();

            eukHess_.template block<vectorPotDim, directorDim>(indexJ, indexI)
                = eukHess_.template block<directorDim, vectorPotDim>(indexI, indexJ).transpose();
          }
        }
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
      Dune::BlockVector<typename std::remove_cvref_t<decltype(mNodal[0])>::template Rebind<ScalarType>::type> mN(fe0.size());
      Dune::BlockVector<typename std::remove_cvref_t<decltype(ANodal[0])>::template Rebind<ScalarType>::type> AN(fe1.size());

      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN[i]=mNodal[globalIndex[1]];
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN[i] = ANodal[globalIndex[1]];
      }

      const auto geo = localView_.element().geometry();

      const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
      Ikarus::ProjectionBasedLocalFunction magnetLocalFunction(localBasisMag, mN);
      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisVecPot, AN);
      using namespace Ikarus::DerivativeDirections;
      ScalarType energy = 0.0;
      for (const auto& gp : rule) {
        localBasisMag.template evaluateFunctionAndJacobian(gp.position(), Nm, dNm);
        localBasisVecPot.template evaluateFunctionAndJacobian(gp.position(), NA, dNA);

        const auto J = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        const auto Jinv             = J.inverse().eval();
        const auto normalizedMag = magnetLocalFunction.evaluateFunction(gp.position()).getValue();
        const auto vectorPot = vectorPotLocalFunction.evaluateFunction(gp.position()).getValue();

        const Eigen::Matrix<double, directorDim, Traits::mydim> gradm
            = magnetLocalFunction.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv));
        const Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
            = vectorPotLocalFunction.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv));

        const Eigen::Vector<ScalarType, 3> curlA(gradA(0, 1), -gradA(0, 0), 0);
        const ScalarType divA = 0;

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);
        energy += (0.5 * (gradm.transpose() * gradm).trace() - 2 * normalizedMag.dot(Hbar) / material.ms)
                  * geo.integrationElement(gp.position()) * gp.weight();  // exchange and zeeman energy

        energy += (0.5 * curlA.squaredNorm() - normalizedMag.dot(curlA) + divA * divA)
                  * geo.integrationElement(gp.position()) * gp.weight();  // demag energy
      }
      return energy;
    }

    mutable Eigen::MatrixXd hEuk;
    mutable Eigen::VectorXd eukGrad;
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
    unsigned int order{};
  };

}  // namespace Ikarus::FiniteElements

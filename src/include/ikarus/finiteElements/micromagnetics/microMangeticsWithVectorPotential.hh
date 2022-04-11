
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

#include <ikarus/finiteElements/interface/feTraits.hh>
#include <ikarus/finiteElements/interface/finiteElementFunctionConcepts.hh>
#include <ikarus/finiteElements/interface/interfaceFiniteElement.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/localFunctions/projectionBasedLocalFunction.hh>
#include <ikarus/localFunctions/standardLocalFunction.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/utils/linearAlgebraTypedefs.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/manifolds/unitVector.hh>

namespace Ikarus {

  namespace Impl {
    template <typename Derived>
    Eigen::Vector<typename Derived::Scalar, 3> jacobianToCurl(const Eigen::MatrixBase<Derived>& jaco) {
      Eigen::Vector<typename Derived::Scalar, 3> curl;
      constexpr int FieldSize = Derived::RowsAtCompileTime;
      constexpr int worlddim  = Derived::ColsAtCompileTime;
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
    using ResultRequirementsType = ResultRequirements<MultiTypeVector>;
    using LocalViewEmbedded      = typename BasisEmbedded::LocalView;
    using LocalViewReduced       = typename BasisReduced::LocalView;

    template <typename VolumeLoad>
    MicroMagneticsWithVectorPotential(BasisEmbedded& globalBasis, BasisReduced& globalBasisRed,
                                      const typename LocalViewEmbedded::Element& element, MagneticMaterial& p_material,
                                      const VolumeLoad& p_volumeLoad, bool p_isInside = true)
        : localView_{globalBasis.localView()},
          localViewReduced{globalBasisRed.localView()},
          volumeLoad(p_volumeLoad),
          material{p_material},
          isInside{p_isInside} {
      localView_.bind(element);
      localViewReduced.bind(element);
      order         = 2 * localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis().order();
      localBasisMag = Ikarus::LocalBasis(localView_.tree().child(Dune::Indices::_0, 0).finiteElement().localBasis());
      localBasisMag.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),
                         bindDerivatives(0, 1));
      localBasisVecPot = Ikarus::LocalBasis(localView_.tree().child(Dune::Indices::_1, 0).finiteElement().localBasis());
      localBasisVecPot.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),
                            bindDerivatives(0, 1));
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
      hred.setZero();

      const auto& mNodal = par.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_0];
      const auto& ANodal = par.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_1];
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      const auto& child0 = localView_.tree().child(_0, 0);
      const auto& fe0    = child0.finiteElement();
      const auto& child1 = localView_.tree().child(_1, 0);
      const auto& fe1    = child1.finiteElement();

      Dune::BlockVector<std::remove_cvref_t<decltype(mNodal[0])>> mN(fe0.size());
      Dune::BlockVector<std::remove_cvref_t<decltype(ANodal[0])>> AN(fe1.size());
      for (auto i = 0U; i < fe0.size(); ++i) {
        const auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN[i]                  = mNodal[globalIndex[1]];
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        const auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN[i]                  = ANodal[globalIndex[1]];
      }

      const auto geo = localView_.element().geometry();

      using namespace Ikarus::DerivativeDirections;
      Ikarus::ProjectionBasedLocalFunction mLocalF(localBasisMag, mN);
      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisVecPot, AN);

      for (const auto& [gpIndex, gp] : mLocalF.viewOverIntegrationPoints()) {
        const double weight = geo.integrationElement(gp.position()) * gp.weight();
        localBasisVecPot.template evaluateJacobian(gp.position(), dNA);

        const auto J             = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        const auto Jinv          = J.inverse().eval();
        const auto normalizedMag = mLocalF.evaluateFunction(gpIndex).getValue();
        const auto gradm         = mLocalF.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));
        const auto vectorPot     = vectorPotLocalFunction.evaluateFunction(gpIndex).getValue();

        const auto dNAdx = (dNA * J.inverse()).eval();
        Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
            = vectorPotLocalFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));

        const Eigen::Vector<double, directorDim> curlA = Impl::jacobianToCurl(gradA).template segment<directorDim>(0);

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);
        if (isInside) {
          for (size_t i = 0; i < fe0.size(); ++i) {
            const int indexI = i * directorCorrectionDim;
            const auto WI
                = mLocalF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(i));
            for (size_t j = 0; j < fe0.size(); ++j) {
              const int indexJ = j * directorCorrectionDim;
              const auto WJ
                  = mLocalF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(j));

              const auto mddHbar = mLocalF.evaluateDerivative(
                  gpIndex, wrt(coeffs, coeffs), along(2 * Hbar / material.ms + curlA), coeffIndices(i, j));

              Eigen::Matrix<double, directorCorrectionDim, directorCorrectionDim> tmp;
              tmp.setZero();
              for (int k = 0; k < Traits::mydim; ++k)
                tmp += WI[k] * WJ[k].transpose()
                       + mLocalF.evaluateDerivative(gpIndex, wrt(spatial(k), coeffs, coeffs), along(gradm.col(k)),
                                                    transformWith(Jinv), coeffIndices(i, j));

              hred.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ) += tmp * weight;
              hred.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ) -= mddHbar * weight;
            }
          }
        }
        const int magEukSize = fe0.size() * directorCorrectionDim;
        for (size_t i = 0; i < fe1.size(); ++i) {
          const int indexI                                                    = magEukSize + i * vectorPotDim;
          const Eigen::Matrix<double, directorDim, vectorPotDim> gradCurlA_dI = Impl::vecToSkewMatrix(dNAdx.row(i));
          ;
          for (size_t j = 0; j < fe1.size(); ++j) {
            const int indexJ                                                    = magEukSize + j * vectorPotDim;
            const Eigen::Matrix<double, directorDim, vectorPotDim> gradCurlA_dJ = Impl::vecToSkewMatrix(dNAdx.row(j));
            ;
            hred.template block<vectorPotDim, vectorPotDim>(indexI, indexJ)
                += (gradCurlA_dI.transpose() * gradCurlA_dJ) * geo.integrationElement(gp.position()) * gp.weight();
          }
        }

        for (size_t i = 0; i < fe0.size(); ++i) {
          const int indexI = i * directorCorrectionDim;
          const auto PmI   = mLocalF.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
          for (size_t j = 0; j < fe1.size(); ++j) {
            const int indexJ                                                    = magEukSize + j * vectorPotDim;
            const Eigen::Matrix<double, directorDim, vectorPotDim> gradCurlA_dJ = Impl::vecToSkewMatrix(dNAdx.row(j));
            if (isInside)
              hred.template block<directorCorrectionDim, vectorPotDim>(indexI, indexJ)
                  -= (PmI * gradCurlA_dJ) * geo.integrationElement(gp.position()) * gp.weight();

            hred.template block<vectorPotDim, directorCorrectionDim>(indexJ, indexI)
                = hred.template block<directorCorrectionDim, vectorPotDim>(indexI, indexJ).transpose();
          }
        }
      }
    }

    [[nodiscard]] int size() const { return localViewReduced.size(); }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& rieGrad) const {
      rieGrad.setZero();
      using namespace Dune::Indices;
      const auto& mNodal = par.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_0];
      const auto& ANodal = par.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_1];
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

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
        AN[i]            = ANodal[globalIndex[1]];
      }

      const auto geo = localView_.element().geometry();
      using namespace Ikarus::DerivativeDirections;
      Ikarus::ProjectionBasedLocalFunction magnetLocalFunction(localBasisMag, mN);
      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisVecPot, AN);

      for (const auto& [gpIndex, gp] : magnetLocalFunction.viewOverIntegrationPoints()) {
        const double weight = geo.integrationElement(gp.position()) * gp.weight();
        localBasisVecPot.evaluateJacobian(gp.position(), dNA);
        const auto J             = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        const auto Jinv          = J.inverse().eval();
        const auto normalizedMag = magnetLocalFunction.evaluateFunction(gpIndex).getValue();
        const auto vectorPot     = vectorPotLocalFunction.evaluateFunction(gpIndex).getValue();

        const auto dNAdx = (dNA * Jinv).eval();
        const Eigen::Matrix<double, directorDim, Traits::mydim> gradm
            = magnetLocalFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));
        const Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
            = vectorPotLocalFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));

        const Eigen::Vector<double, directorDim> curlA = Impl::jacobianToCurl(gradA).template segment<directorDim>(0);

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);

        if (isInside) {
          for (size_t i = 0; i < fe0.size(); ++i) {
            const int index = i * directorCorrectionDim;
            const auto WI   = magnetLocalFunction.evaluateDerivative(gpIndex, wrt(spatialall, coeffs),
                                                                   transformWith(Jinv), coeffIndices(i));
            const auto PmI  = magnetLocalFunction.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
            Eigen::Vector<double, directorCorrectionDim> tmp;
            tmp.setZero();
            for (int j = 0; j < Traits ::mydim; ++j)
              tmp += WI[j] * gradm.col(j);

            rieGrad.template segment<directorCorrectionDim>(index) += (tmp - PmI * curlA) * weight;
            rieGrad.template segment<directorCorrectionDim>(index) -= (2 * PmI * Hbar / material.ms) * weight;
          }
        }
        const int magEukSize = fe0.size() * directorCorrectionDim;
        for (size_t i = 0; i < fe1.size(); ++i) {
          const int index                                                     = magEukSize + i * vectorPotDim;
          const Eigen::Matrix<double, directorDim, vectorPotDim> gradCurlA_dI = Impl::vecToSkewMatrix(dNAdx.row(i));
          rieGrad.template segment<vectorPotDim>(index)
              += (-gradCurlA_dI.transpose() * (normalizedMag * isInside) + gradCurlA_dI.transpose() * (curlA)) * weight;
        }
      }
    }

    [[nodiscard]] typename Traits::ScalarType calculateScalar(const FERequirementType& par) const {
      Eigen::VectorXd dx(localView_.size());
      dx.setZero();

      return this->calculateScalarImpl(par, dx);
    }

    void calculateAt(const ResultRequirementsType& req, const Eigen::Vector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::Indices;
      const auto& mNodal = req.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_0];
      const auto& ANodal = req.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_1];
      const auto& lambda = req.getParameter(Ikarus::FEParameter::loadfactor);
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
        AN[i]            = ANodal[globalIndex[1]];
      }

      const auto geo = localView_.element().geometry();

      auto gp = toFieldVector(local);
      //      localBasisMag.template evaluateFunctionAndJacobian(gp, Nm, dNm);
      //      localBasisVecPot.template evaluateFunctionAndJacobian(gp, NA, dNA);

      using namespace Ikarus::DerivativeDirections;
      Ikarus::ProjectionBasedLocalFunction magnetLocalFunction(localBasisMag, mN);

      const auto J    = toEigenMatrix(geo.jacobianTransposed(gp)).transpose().eval();
      const auto Jinv = J.inverse().eval();
      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisVecPot, AN);

      const auto normalizedMag = magnetLocalFunction.evaluateFunction(gp).getValue();
      const auto vectorPot     = vectorPotLocalFunction.evaluateFunction(gp).getValue();

      const Eigen::Matrix<double, directorDim, Traits::mydim> gradm
          = magnetLocalFunction.evaluateDerivative(gp, wrt(spatialall), transformWith(Jinv));
      const Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
          = vectorPotLocalFunction.evaluateDerivative(gp, wrt(spatialall), transformWith(Jinv));

      const Eigen::Vector<double, 3> curlA          = Impl::jacobianToCurl(gradA);
      const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp), lambda);

      typename ResultTypeMap<double>::ResultArray resv;
      if( req.isResultRequested( ResultType::gradientNormOfMagnetization)) {
        resv.resize(1,1);
        resv(0,0)=isInside ? gradm.norm() : 0.0;
        result.insertOrAssignResult(ResultType::gradientNormOfMagnetization,resv);
      }
      if( req.isResultRequested( ResultType::BField)) {
        resv = curlA/(material.mu0*material.ms)*std::sqrt(2.0);
        result.insertOrAssignResult(ResultType::BField,resv);
      }
      if( req.isResultRequested( ResultType::HField)) {
        resv.setZero(3,1);
        resv(0,0) = curlA[0]/(Dune::power(material.mu0, 2)*material.ms)*std::sqrt(2.0)
            - normalizedMag[0]*material.ms*static_cast<double>(isInside);
        resv(1,0) = curlA[1]/(Dune::power(material.mu0, 2)*material.ms)*std::sqrt(2.0)
            - normalizedMag[1]*material.ms*static_cast<double>(isInside);
        if constexpr (directorDim==3)
          resv(2,0) = curlA[2]/(Dune::power(material.mu0, 2)*material.ms)*std::sqrt(2.0)
              - normalizedMag[2]*material.ms*static_cast<double>(isInside);
        result.insertOrAssignResult(ResultType::HField,resv);
      }

    }

  private:
    template <class ScalarType>
    ScalarType calculateScalarImpl(const FERequirementType& par, Eigen::VectorX<ScalarType>& dx_) const {
      using namespace Dune::Indices;
      const auto& mNodal = par.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_0];
      const auto& ANodal = par.getSolution(Ikarus::FESolutions::magnetizationAndVectorPotential)[_1];
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      auto& child0    = localView_.tree().child(_0, 0);
      const auto& fe0 = child0.finiteElement();
      auto& child1    = localView_.tree().child(_1, 0);

      const auto& fe1 = child1.finiteElement();
      Dune::BlockVector<typename std::remove_cvref_t<decltype(mNodal[0])>::template Rebind<ScalarType>::type> mN(
          fe0.size());
      Dune::BlockVector<typename std::remove_cvref_t<decltype(ANodal[0])>::template Rebind<ScalarType>::type> AN(
          fe1.size());

      for (auto i = 0U; i < fe0.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_0, 0).localIndex(i));
        mN[i]            = mNodal[globalIndex[1]];
      }

      const int magElementEntries = fe0.size() * directorDim;
      for (auto i = 0U; i < fe1.size(); ++i) {
        auto globalIndex = localView_.index(localView_.tree().child(_1, 0).localIndex(i));
        AN[i]            = ANodal[globalIndex[1]];
      }

      const auto geo = localView_.element().geometry();

      const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
      Ikarus::ProjectionBasedLocalFunction magnetLocalFunction(localBasisMag, mN);
      Ikarus::StandardLocalFunction vectorPotLocalFunction(localBasisVecPot, AN);
      using namespace Ikarus::DerivativeDirections;
      ScalarType energy = 0.0;
      for (const auto& gp : rule) {
        const auto J             = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
        const auto Jinv          = J.inverse().eval();
        const auto normalizedMag = magnetLocalFunction.evaluateFunction(gp.position()).getValue();
        const auto vectorPot     = vectorPotLocalFunction.evaluateFunction(gp.position()).getValue();

        const Eigen::Matrix<double, directorDim, Traits::mydim> gradm
            = magnetLocalFunction.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv));
        const Eigen::Matrix<double, vectorPotDim, Traits::mydim> gradA
            = vectorPotLocalFunction.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv));

        const Eigen::Vector<ScalarType, directorDim> curlA
            = Impl::jacobianToCurl(gradA).template segment<directorDim>(0);

        const Eigen::Vector<double, directorDim> Hbar = volumeLoad(toEigenVector(gp.position()), lambda);
        if (isInside)
          energy += (0.5 * (gradm.transpose() * gradm).trace() - 2 * normalizedMag.dot(Hbar) / material.ms + 0.5 * 1)
                    * geo.integrationElement(gp.position()) * gp.weight();  // exchange and zeeman energy

        energy += (0.5 * curlA.squaredNorm() - isInside * normalizedMag.dot(curlA))
                  * geo.integrationElement(gp.position()) * gp.weight();  // demag energy
      }
      return energy;
    }

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
    bool isInside{};
  };

}  // namespace Ikarus::FiniteElements

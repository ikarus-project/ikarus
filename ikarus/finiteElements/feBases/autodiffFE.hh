// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

////
//
#pragma once
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {
  template <typename RealElement,
            bool forceAutoDiff = false>
  class AutoDiffFE : public RealElement {
  public:
    using Base              = RealElement;
    using Basis             = Base::Basis;
    using LocalView         = typename Basis::FlatBasis::LocalView;
    using Traits            = TraitsFromLocalView<LocalView, RealElement::useEigenRef>;
    using Element           = typename LocalView::Element;
    using FERequirementType = typename RealElement::FERequirementType;

    void fillReducedH(const auto& par,const auto& h,const auto& g, auto& hred)const
    {
//      std::cout<<"h\n"<<std::endl;
//      std::cout<<h<<std::endl;
//      std::cout<<"g\n"<<std::endl;
//      std::cout<<g.transpose()<<std::endl;
      if constexpr (not std::is_same_v<FERequirementType,FErequirements<>>)
      {
        using namespace Dune::Indices;
        const auto &mNodal = par.getGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector)[_0];
        const auto &dNodal = par.getGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector)[_1];
        const auto &child0 = this->localViewBlocked_.tree().child(_0, 0);
        const auto &child1 = this->localViewBlocked_.tree().child(_1, 0);
        const auto &fe0    = child0.finiteElement();
        const auto &fe1    = child1.finiteElement();
        Dune::BlockVector<std::remove_cvref_t<decltype(mNodal[0])>> localMidSurfaceConfiguration(fe0.size());
        for (auto i = 0U; i < fe0.size(); ++i) {
          const auto globalIndex = this->localViewBlocked().index(this->localViewBlocked().tree().child(_0, 0).localIndex(i));
          localMidSurfaceConfiguration[i] = mNodal[globalIndex[1]];
        }

        Dune::BlockVector<std::remove_cvref_t<decltype(dNodal[0])>> localDirectorConfiguration(fe1.size());
        for (auto i = 0U; i < fe1.size(); ++i) {
          const auto globalIndex = this->localViewBlocked().index(this->localViewBlocked().tree().child(_1, 0).localIndex(i));
          localDirectorConfiguration[i] = dNodal[globalIndex[1]];
        }
        using Manifold0 = std::remove_cvref_t<decltype(localMidSurfaceConfiguration[0])>;
        using Manifold1 = std::remove_cvref_t<decltype(localDirectorConfiguration[0])>;
        constexpr int valueSize0 =Manifold0::valueSize;
        constexpr int valueSize1 =Manifold1::valueSize;
        constexpr int correctionSize0  =Manifold0::correctionSize ;
        constexpr int correctionSize1  =Manifold1::correctionSize ;
        const int nDofs0 = this->localViewBlocked().tree().child(_0, 0).finiteElement().size()*correctionSize0;
        const int nDofs1 = this->localViewBlocked().tree().child(_0, 0).finiteElement().size()*correctionSize1;
        hred.block(0,0,nDofs0,nDofs0)= h.block(0,0,nDofs0,nDofs0);
//        std::cout<<"hredBB\n"<<std::endl;
//        std::cout<<hred<<std::endl;

        for (int i = 0; i < localDirectorConfiguration.size(); i++) {

          const int index3i = nDofs0+valueSize1*i;
          const int index2i = nDofs0+correctionSize1*i;
          const auto BLAIT= localDirectorConfiguration[i].orthonormalFrame().transpose().eval();
          for (int j = 0; j < localDirectorConfiguration.size(); j++) {
            const int index3j = nDofs0+valueSize1*j;
            const int index2j = nDofs0+correctionSize1*j;
            const auto BLAIJ= localDirectorConfiguration[j].orthonormalFrame();
            hred.template block<correctionSize1,correctionSize1>(index2i,index2j) =BLAIT*h.template block<valueSize1,valueSize1>(index3i,index3j)*BLAIJ;

          }
          for (int j = 0; j < localMidSurfaceConfiguration.size(); j++) {
            const int index3j = valueSize0*j;
            const int index2j = correctionSize0*j;
            hred.template block<correctionSize1,correctionSize0>(index2i,index2j) =BLAIT*h.template block<valueSize1,valueSize0>(index3i,index3j);

          }
            hred.template block<correctionSize1,correctionSize1>(index2i,index2i)+=localDirectorConfiguration[i].weingarten(g.template segment<valueSize1>(index3i));
        }
//        std::cout<<"hredB\n"<<std::endl;
//        std::cout<<hred<<std::endl;
        for (int i = 0; i < localMidSurfaceConfiguration.size(); i++) {
          const int index3i = valueSize0*i;
          const int index2i = correctionSize0*i;
          for (int j = 0; j < localDirectorConfiguration.size(); j++) {
            const int index3j = nDofs0+valueSize1*j;
            const int index2j = nDofs0+correctionSize1*j;
            const auto BLAIJ= localDirectorConfiguration[j].orthonormalFrame();
            hred.template block<correctionSize0,correctionSize1>(index2i,index2j) =h.template block<valueSize0,valueSize1>(index3i,index3j)*BLAIJ;

          }
          }
         hred.block(0,nDofs0,nDofs0,nDofs1)=hred.block(nDofs0,0,nDofs1,nDofs0).transpose();
        std::cout<<"hredA\n"<<std::endl;
        std::cout<<hred<<std::endl;
      }else
        hred=h;
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> hred) const {
      if constexpr (requires { RealElement::calculateMatrix(par, hred); } and not forceAutoDiff) {
        RealElement::calculateMatrix(par, hred);
      } else if constexpr (requires {
                             this->template calculateVectorImpl<autodiff::dual>(
                                 par, std::declval<typename Traits::template VectorType<autodiff::dual>>(),
                                 std::declval<const Eigen::VectorXdual&>());
                           }) {
        /// This is only valid if the external forces are independent of displacements, for e.g., no follower forces are
        /// applied
        Eigen::VectorXdual dx(this->ndofEmbedded());
        Eigen::VectorXdual g(this->ndofEmbedded());
        Eigen::MatrixXd h(this->ndofEmbedded(),this->ndofEmbedded());
        dx.setZero();
        auto f = [&](auto& x) -> auto& {
          g.setZero();
          this->template calculateVectorImpl<autodiff::dual>(par, g, x);
          return g;
        };
        jacobian(f, autodiff::wrt(dx), at(dx), g, h);
        fillReducedH(par,h,g,hred);
      } else if constexpr (requires {
                             this->template calculateScalarImpl<autodiff::dual2nd>(
                                 par, std::declval<typename Traits::template VectorType<autodiff::dual2nd>>());
                           }) {
        Eigen::VectorXdual2nd dx(this->ndofEmbedded());
        Eigen::VectorXd g;
        Eigen::MatrixXd h(this->ndofEmbedded(),this->ndofEmbedded());
        autodiff::dual2nd e;
        dx.setZero();
        auto f = [&](auto& x) { return this->template calculateScalarImpl<autodiff::dual2nd>(par, x); };
        hessian(f, autodiff::wrt(dx), at(dx), e, g, h);
        fillReducedH(par,h,g,hred);
      } else
        static_assert(Ikarus::Std::DummyFalse<AutoDiffFE>::value,
                      "Appropriate calculateScalarImpl or calculateVectorImpl functions are not implemented for the "
                      "chosen element.");
    }

    inline void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> gRed) const {
      if constexpr (requires {
                      this->template calculateVectorImpl<double>(
                          par, std::declval<typename Traits::template VectorType<double>>(),
                          std::declval<const Eigen::VectorXd&>());
                    }
                    and not forceAutoDiff) {
        return this->template calculateVectorImpl<double>(par, gRed);
      } else if constexpr (requires {
                             this->template calculateScalarImpl<autodiff::dual>(
                                 par, std::declval<const Eigen::VectorXdual&>());
                           }) {
        Eigen::VectorXdual dx(this->ndofEmbedded());
        Eigen::VectorXd g(this->ndofEmbedded());
//        Eigen::VectorXdual gRed(this->size());
        dx.setZero();
        gRed.setZero();
        autodiff::dual e;
        auto f = [&](auto& x) { return this->template calculateScalarImpl<autodiff::dual>(par, x); };
        gradient(f, autodiff::wrt(dx), at(dx), e, g);
        if constexpr (not std::is_same_v<FERequirementType,FErequirements<>>)
        {
          using namespace Dune::Indices;
          const auto &dNodal = par.getGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector)[_1];
          const auto &child1 = this->localViewBlocked_.tree().child(_1, 0);
          const auto &fe1    = child1.finiteElement();
          Dune::BlockVector<std::remove_cvref_t<decltype(dNodal[0])>> localDirectorConfiguration(fe1.size());
          for (auto i = 0U; i < fe1.size(); ++i) {
            const auto globalIndex = this->localViewBlocked().index(this->localViewBlocked().tree().child(_1, 0).localIndex(i));
              localDirectorConfiguration[i] = dNodal[globalIndex[1]];
          }
          const int nDofs0 = this->localViewBlocked().tree().child(_0, 0).finiteElement().size();
          gRed.segment(0,nDofs0)= g.segment(0,nDofs0);
          for (int i = 0; i < localDirectorConfiguration.size(); i++) {
             const int index3 = nDofs0+3*i;
             const int index2 = nDofs0+2*i;
             gRed.template segment<2>(index2)= localDirectorConfiguration[i].orthonormalFrame().transpose()*g.template segment<3>(index3);
          }
        }else
          gRed=g;
      } else
        static_assert(Ikarus::Std::DummyFalse<AutoDiffFE>::value,
                      "Appropriate calculateScalarImpl function is not implemented for the "
                      "chosen element.");
    }

    void calculateLocalSystem(const FERequirementType& par, typename Traits::template MatrixType<> h,
                              typename Traits::template VectorType<> g) const {
      Eigen::VectorXdual2nd dx(this->ndofEmbedded());
      dx.setZero();
      auto f = [&](auto& x) { return this->calculateScalarImpl(par, x); };
      hessian(f, autodiff::wrt(dx), at(dx), g, h);
    }

    [[nodiscard]] double calculateScalar(const FERequirementType& par) const {
      if constexpr (requires { RealElement::calculateScalar(par); }) {
        return RealElement::calculateScalar(par);
      } else if constexpr (requires { this->calculateScalarImpl(par); }) {
        return this->calculateScalarImpl(par);
      } else
        static_assert(Ikarus::Std::DummyFalse<AutoDiffFE>::value,
                      "Appropriate calculateScalar and calculateScalarImpl functions are not implemented for the "
                      "chosen element.");
    }

    const RealElement& getFE() const { return *this; }

    template <typename... Args>
    explicit AutoDiffFE(Args&&... args) : RealElement{std::forward<Args>(args)...} {}
  };
}  // namespace Ikarus


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

#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/localFunctions/meta.hh>

namespace Ikarus{

  static constexpr int maxEASParameter3d = 21;
  static constexpr int maxEASParameter2d = 8;

  enum class EASType{
    none,
    Q1E4,
    Q1E5,
    Q1E7,
    H1E9,
    H1E21
  };

  template<typename Geometry>
  Eigen::Matrix3d calcTransformationMatrix2D(const Geometry& geometry){

    const auto quadPos0 = geometry.type().template geometry<0>(0).center();

    const auto jacobianinvT0 = toEigenMatrix(geometry.jacobianInverseTransposed(quadPos0)); //J^{-1}.Transpose() in Dune = J^{-1}
    const auto detJ0 = geometry.integrationElement(quadPos0); //determinant(J)

    auto jaco = (jacobianinvT0).inverse().eval();
    auto J11 = jaco(0,0);
    auto J12 = jaco(0,1);
    auto J21 = jaco(1,0);
    auto J22 = jaco(1,1);

    Eigen::Matrix3d T0;
    T0 << J11*J11 , J12*J12 , J11*J12 ,
        J21*J21 , J22*J22 , J21*J22 ,
        2.0*J11*J21 , 2.0*J12*J22 , J21*J12 + J11*J22;

    return T0.inverse() * detJ0;
  }

    template<typename Geometry>
    Eigen::Matrix3d calcTransformationMatrix3D(const Geometry& geometry){

        const auto quadPos0 = geometry.type().template geometry<0>(0).center();

        const auto jacobianinvT0 = toEigenMatrix(geometry.jacobianInverseTransposed(quadPos0)); //J^{-1}.Transpose() in Dune = J^{-1}
        const auto detJ0 = geometry.integrationElement(quadPos0); //determinant(J)

        auto jaco = (jacobianinvT0).inverse().eval();
        auto J11 = jaco(0,0);
        auto J12 = jaco(0,1);
        auto J13 = jaco(0,2);
        auto J21 = jaco(1,0);
        auto J22 = jaco(1,1);
        auto J23 = jaco(1,2);
        auto J31 = jaco(2,0);
        auto J32 = jaco(2,1);
        auto J33 = jaco(2,2);

        Eigen::Matrix<double,6,6> T0;
        T0 << J11*J11 , J12*J12 , J13*J13 , J11*J12 , J11*J13 , J12*J13 ,
                J21*J21 , J22*J22 , J23*J23 , J21*J22 , J21*J23 , J22*J23 ,
                J31*J31 , J32*J32 , J33*J33 , J31*J32 , J31*J33 , J32*J33 ,
                2.0*J11*J21 , 2.0*J12*J22 , 2.0*J13*J23 , J11*J22 + J21*J12 , J11*J23 + J21*J13 , J12*J23 + J22*J13 ,
                2.0*J11*J31 , 2.0*J12*J32 , 2.0*J13*J33 , J11*J32 + J31*J12 , J11*J33 + J31*J13 , J12*J33 + J32*J13 ,
                2.0*J31*J21 , 2.0*J32*J22 , 2.0*J33*J23 , J31*J22 + J21*J32 , J31*J23 + J21*J33 , J32*J23 + J22*J33 ;

        return T0.inverse() * detJ0;
    }


  template<typename Geometry>
  struct EASQ1E4
  {
    static constexpr int strainSize = 3;
    static constexpr int enhancedStrainSize = 4;

    EASQ1E4() = default;
    EASQ1E4(const Geometry& geometry)
    : geometry{std::make_unique<Geometry>(geometry)},
          T0InverseTransformed{calcTransformationMatrix2D(geometry)}
    {}

    auto calcM(const Dune::FieldVector<double,2>& quadPos) const
    {
      Eigen::Matrix<double,strainSize,enhancedStrainSize> M;
      M.setZero(strainSize,enhancedStrainSize);
      const double xi = quadPos[0];
      const double eta = quadPos[1];
      M(0,0) = xi-0.5;
      M(1,1) = eta-0.5;
      M(2,2) = xi-0.5;
      M(2,3) = eta-0.5;
      const double detJ = geometry->integrationElement(quadPos);
      M = T0InverseTransformed/detJ * M ;
      return M;
    }

     std::unique_ptr<Geometry> geometry;
    Eigen::Matrix3d T0InverseTransformed;
  };

  template<typename Geometry>
  struct EASQ1E5
  {
    static constexpr int strainSize = 3;
    static constexpr int enhancedStrainSize = 5;

    EASQ1E5() = default;
    EASQ1E5(const Geometry& geometry)
        : geometry{std::make_unique<Geometry>(geometry)}
          ,T0InverseTransformed{calcTransformationMatrix2D(geometry)}
    {}

    auto calcM(const Dune::FieldVector<double,2>& quadPos) const
    {
      Eigen::Matrix<double,strainSize,enhancedStrainSize> M;
      M.setZero();
      const double xi = quadPos[0];
      const double eta = quadPos[1];
      M(0,0) = xi-0.5;
      M(1,1) = eta-0.5;
      M(2,2) = xi-0.5;
      M(2,3) = eta-0.5;
      M(2,4) = (xi-0.5)*(eta-0.5);
      const double detJ = geometry->integrationElement(quadPos);
      M = T0InverseTransformed/detJ * M ;
      return M;
    }

    std::unique_ptr<Geometry> geometry;
    Eigen::Matrix3d T0InverseTransformed;
  };

    template<typename Geometry>
    struct EASQ1E7
    {
        static constexpr int strainSize = 3;
        static constexpr int enhancedStrainSize = 7;

        EASQ1E7() = default;
        EASQ1E7(const Geometry& geometry)
                : geometry{std::make_unique<Geometry>(geometry)}
                ,T0InverseTransformed{calcTransformationMatrix2D(geometry)}
        {}

        auto calcM(const Dune::FieldVector<double,2>& quadPos) const
        {
            Eigen::Matrix<double,strainSize,enhancedStrainSize> M;
            M.setZero();
            const double xi = quadPos[0];
            const double eta = quadPos[1];
            M(0,0) = xi-0.5;
            M(1,1) = eta-0.5;
            M(2,2) = xi-0.5;
            M(2,3) = eta-0.5;
            M(0,4) = (xi-0.5)*(eta-0.5);
            M(1,5) = (xi-0.5)*(eta-0.5);
            M(2,6) = (xi-0.5)*(eta-0.5);
            const double detJ = geometry->integrationElement(quadPos);
            M = T0InverseTransformed/detJ * M ;
            return M;
        }

        std::unique_ptr<Geometry> geometry;
        Eigen::Matrix3d T0InverseTransformed;
    };

    template<typename Geometry>
    struct EASH1E9
    {
        static constexpr int strainSize = 6;
        static constexpr int enhancedStrainSize = 9;

        EASH1E9() = default;
        EASH1E9(const Geometry& geometry)
                : geometry{std::make_unique<Geometry>(geometry)}
                ,T0InverseTransformed{calcTransformationMatrix3D(geometry)}
        {}

        auto calcM(const Dune::FieldVector<double,3>& quadPos) const
        {
            Eigen::Matrix<double,strainSize,enhancedStrainSize> M;
            M.setZero();
            const double xi = quadPos[0];
            const double eta = quadPos[1];
            const double zeta = quadPos[2];
            M(0,0) = xi-0.5;
            M(1,1) = eta-0.5;
            M(2,2) = zeta-0.5;
            M(3,3) = xi-0.5;
            M(3,4) = eta-0.5;
            M(4,5) = xi-0.5;
            M(4,6) = zeta-0.5;
            M(5,7) = eta-0.5;
            M(5,8) = zeta-0.5;
            const double detJ = geometry->integrationElement(quadPos);
            M = T0InverseTransformed/detJ * M ;
            return M;
        }
        std::unique_ptr<Geometry> geometry;
        Eigen::Matrix<double,6,6> T0InverseTransformed;
    };

    template<typename Geometry>
    struct EASH1E21
    {
        static constexpr int strainSize = 6;
        static constexpr int enhancedStrainSize = 21;

        EASH1E21() = default;
        EASH1E21(const Geometry& geometry)
                : geometry{std::make_unique<Geometry>(geometry)}
                ,T0InverseTransformed{calcTransformationMatrix3D(geometry)}
        {}

        auto calcM(const Dune::FieldVector<double,3>& quadPos) const
        {
            Eigen::Matrix<double,strainSize,enhancedStrainSize> M;
            M.setZero();
            const double xi = quadPos[0];
            const double eta = quadPos[1];
            const double zeta = quadPos[2];
            M(0,0) = xi-0.5;
            M(1,1) = eta-0.5;
            M(2,2) = zeta-0.5;
            M(3,3) = xi-0.5;
            M(3,4) = eta-0.5;
            M(4,5) = xi-0.5;
            M(4,6) = zeta-0.5;
            M(5,7) = eta-0.5;
            M(5,8) = zeta-0.5;

            M(3,9) = (xi-0.5)*(zeta-0.5);
            M(3,10) = (eta-0.5)*(zeta-0.5);
            M(4,11) = (xi-0.5)*(eta-0.5);
            M(4,12) = (eta-0.5)*(zeta-0.5);
            M(5,13) = (xi-0.5)*(eta-0.5);
            M(5,14) = (xi-0.5)*(zeta-0.5);

            M(0,15) = (xi-0.5)*(eta-0.5);
            M(0,16) = (xi-0.5)*(zeta-0.5);
            M(1,17) = (xi-0.5)*(eta-0.5);
            M(1,18) = (eta-0.5)*(zeta-0.5);
            M(2,19) = (xi-0.5)*(zeta-0.5);
            M(2,20) = (eta-0.5)*(zeta-0.5);

            const double detJ = geometry->integrationElement(quadPos);
            M = T0InverseTransformed/detJ * M ;
            return M;
        }

        std::unique_ptr<Geometry> geometry;
        Eigen::Matrix<double,6,6> T0InverseTransformed;
    };

  /// 2D - Q1 (4-node) , 3D - H1 (8-node) variants
  template<typename Geometry>
  using EAS2DVariant = std::variant<EASQ1E4<Geometry>,EASQ1E5<Geometry>,EASQ1E7<Geometry>>;
  template<typename Geometry>
  using EAS3DVariant = std::variant<EASH1E9<Geometry>,EASH1E21<Geometry>>;

  template<typename DisplacementBasedElement>
class EnhancedAssumedStrains : public DisplacementBasedElement {
  public:

    using FERequirementType = FErequirements<Eigen::VectorXd>;
    using LocalView         = typename DisplacementBasedElement::LocalView;
    using GridView         = typename DisplacementBasedElement::GridView;
    using Traits= typename DisplacementBasedElement::Traits;
    using GlobalIndex= typename DisplacementBasedElement::GlobalIndex;

    template <typename Basis,typename VolumeLoad, typename NeumannBoundaryLoad>
    EnhancedAssumedStrains(Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu, const BoundaryPatch<GridView> * neumannBoundary,
                           const NeumannBoundaryLoad& neumannBoundaryLoad,
                           const VolumeLoad& p_volumeLoad)
        :DisplacementBasedElement(globalBasis,element,emod,nu,neumannBoundary,neumannBoundaryLoad,p_volumeLoad)
    {

    }

    double calculateScalar(const FERequirementType& par)const {

      DUNE_THROW(Dune::NotImplemented, "EAS element do not support any scalar calculations");
      return 0.0;
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      DisplacementBasedElement::calculateVector(par,g);
}

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      using namespace DerivativeDirections;

        /// fill h with displacement-based stiffness.
        /// It is assumed to be assembled block-wise on element level.
        /// This means the displacements x,y,z of node I are grouped together.

      DisplacementBasedElement::calculateMatrix(par,h);

      if(onlyDisplacementBase)
        return ;

      auto strainFunction = DisplacementBasedElement::getStrainFunction(par);
      auto C = DisplacementBasedElement::getMaterialTangentFunction(par);
      auto& localView = DisplacementBasedElement::getLocalView();
      auto geo = localView.element().geometry();

      auto& first_child = localView.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      assert(((fe.size()== 4 and Traits::mydim==2) or (fe.size()== 8 and Traits::mydim==3)) && "EAS only supported for Q1 or H1 elements");
      L.setZero(enhancedStrainSize,localView.size());
      D.setZero();
      for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
        std::visit([&](const auto& easfunction){ M = easfunction.calcM(gp.position()); },easVariant);
        const auto Jinv = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().inverse().eval();
        const auto Ceval = C(gpIndex);
        const double detJ = geo.integrationElement(gp.position());
        D += M.transpose() * Ceval * M * detJ * gp.weight();

      for (size_t i = 0; i < fe.size(); ++i) {
        const size_t I    = Traits::worlddim * i;
        const auto Bi = strainFunction.evaluateDerivative(gpIndex,wrt( coeff(i) ), transformWith(Jinv));

        L.block(0,I,enhancedStrainSize,Traits::worlddim) += M.transpose() * Ceval * Bi * detJ * gp.weight();

        }
      }
      h-= L.transpose()*D.inverse()*L; //exploit symmetry
    }

void setEASType(int numberOfEASParameters)
    {
  enhancedStrainSize=numberOfEASParameters;
  D.setZero(enhancedStrainSize,enhancedStrainSize);

  if constexpr (Traits::mydim == 2) {
      switch (numberOfEASParameters) {
          case 0:
              onlyDisplacementBase = true;
              break;
          case 4:
              easVariant = EASQ1E4(DisplacementBasedElement::getLocalView().element().geometry());
              break;
          case 5:
              easVariant = EASQ1E5(DisplacementBasedElement::getLocalView().element().geometry());
              break;
          case 7:
              easVariant = EASQ1E7(DisplacementBasedElement::getLocalView().element().geometry());
              break;
          default:
              DUNE_THROW(Dune::NotImplemented, "The given EAS parameters are not available.");
              break;
      }
  }
    else if constexpr (Traits::mydim == 3) {
        switch (numberOfEASParameters) {
            case 0:
                onlyDisplacementBase = true;
                break;
            case 9:
                easVariant = EASH1E9(DisplacementBasedElement::getLocalView().element().geometry());
                break;
            case 21:
                easVariant = EASH1E21(DisplacementBasedElement::getLocalView().element().geometry());
                break;
            default:
                DUNE_THROW(Dune::NotImplemented, "The given EAS parameters are not available.");
                break;
        }
    }
}

private:
  using EAS2DVariantImpl = EAS2DVariant<typename LocalView::Element::Geometry>;
  using EAS3DVariantImpl = EAS3DVariant<typename LocalView::Element::Geometry>;
  std::conditional_t<Traits::mydim == 2,EAS2DVariantImpl,EAS3DVariantImpl> easVariant;
  mutable Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M;
  mutable Eigen::MatrixXd D;
  mutable Eigen::MatrixXd L;
  int enhancedStrainSize;
  bool onlyDisplacementBase{false};
};
}




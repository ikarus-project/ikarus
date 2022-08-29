
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
  Eigen::Matrix3d calcTransformationMatrix2D(const Geometry& geometry)
  {
    Dune::FieldVector<double,2> quadPos0;
    quadPos0[0] = 0.5; // Center of the Element in Domain [0,1]
    quadPos0[1] = 0.5; // Center of the Element in Domain [0,1]

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
  using EAS2dVariant = std::variant<EASQ1E4<Geometry>,EASQ1E5<Geometry>>;
//  using EAS3dVariant = std::variant<EASH1E9,EASH1E21>;


  template<typename DisplacementBasedElement>
class EnhancedAssumedStrains : public DisplacementBasedElement {
  public:



    static constexpr int strainSize = 3;

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
      DisplacementBasedElement::calculateMatrix(par,h); // fill h with displacement-based stiffnesses
      //is assumed to be assembled block-wise on element level. This means the displacements x,y,z of node I are grouped together

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

  switch (numberOfEASParameters) {
    case 0:
      onlyDisplacementBase=true;
      break;
    case 4:
      easVariant= EASQ1E4(DisplacementBasedElement::getLocalView().element().geometry());
      break;
    case 5:
      easVariant= EASQ1E5(DisplacementBasedElement::getLocalView().element().geometry());
      break;
    default:
      DUNE_THROW(Dune::NotImplemented,"The given EAS parameters are not available.");
      break;
  }


}

private:
  EAS2dVariant<typename LocalView::Element::Geometry> easVariant;
  mutable Eigen::Matrix<double,strainSize,Eigen::Dynamic> M;
  mutable Eigen::MatrixXd D;
  mutable Eigen::MatrixXd L;
  int enhancedStrainSize;
  bool onlyDisplacementBase{false};
};




}
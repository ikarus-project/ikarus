
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

namespace Ikarus{


  enum class EASType{
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
  struct EASE4
  {
    static constexpr int strainSize = 3;
    static constexpr int enhancedStrainSize = 4;

    EASE4(const Geometry& geometry)
    : geometry{geometry},
          T0InverseTransformed{calcTransformationMatrix2D(geometry)}
    {}

    auto calcM(const Dune::FieldVector<double,2>& quadPos)
    {
      Eigen::Matrix<double,strainSize,enhancedStrainSize> M;
      M.setZero();
      M(0,0) = quadPos[0]-0.5;
      M(1,1) = quadPos[1]-0.5;
      M(2,2) = quadPos[0]-0.5;
      M(2,3) = quadPos[1]-0.5;
      const double detJ = geometry.integrationElement(quadPos);
      M = T0InverseTransformed/detJ * M ;
      return M;
    }

    const Geometry& geometry;
    Eigen::Matrix3d T0InverseTransformed;
  };

  template<typename Geometry>
  struct EASE5
  {
    static constexpr int strainSize = 3;
    static constexpr int enhancedStrainSize = 5;

    EASE5(const Geometry& geometry)
        : geometry{geometry}
          ,T0InverseTransformed{calcTransformationMatrix2D(geometry)}
    {}

    auto calcM(const Dune::FieldVector<double,2>& quadPos)
    {
      Eigen::Matrix<double,strainSize,enhancedStrainSize> M;
      M.setZero();
      M(0,0) = quadPos[0]-0.5;
      M(1,1) = quadPos[1]-0.5;
      M(2,2) = quadPos[0]-0.5;
      M(2,3) = quadPos[1]-0.5;
      M(2,4) = (quadPos[0]-0.5)*(quadPos[1]-0.5);
      M = T0InverseTransformed  * M ;
      return M;
    }

    const Geometry& geometry;
    Eigen::Matrix3d T0InverseTransformed;
  };



  template<typename DisplacementBasedElement>
class EnhancedAssumedStrains  {
  public:

    static constexpr int strainSize = 3;
    static constexpr int enhancedStrainSize = 4;
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    using LocalView         = typename DisplacementBasedElement::LocalView;
    using GridView         = typename DisplacementBasedElement::GridView;
    using Traits= typename DisplacementBasedElement::Traits;

    template <typename Basis,typename VolumeLoad, typename NeumannBoundaryLoad>
    EnhancedAssumedStrains(Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu, const BoundaryPatch<GridView> * neumannBoundary,
                           const NeumannBoundaryLoad& neumannBoundaryLoad,
                           const VolumeLoad& p_volumeLoad)
        :displacementBasedElement(globalBasis,element,emod,nu,neumannBoundary,neumannBoundaryLoad,p_volumeLoad)
    {

    }

    double calculateScalar(const FERequirementType& par)const {

      DUNE_THROW(Dune::NotImplemented, "EAS element do not support any scalar calculations");
      return 0.0;
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      DUNE_THROW(Dune::NotImplemented, "EAS element do not support any vector calculations. Since it is only used for linear elasticity right now. Later the internal forces need to be calculated here");

}

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      displacementBasedElement.calculateMatrix(par,h); // fill h with displacement-based stiffnesses
      //is assumed to be assembled block-wise on element level. This means the displacements x,y,z of node I are grouped together

      auto strainFunction = displacementBasedElement.getStrainFunction(par);
      auto C = displacementBasedElement.getMaterialTangentFunction(par);
      auto& localView = displacementBasedElement.getLocalView();
      auto geo = localView.element().geometry();
      auto easFunction = EASE4(geo);
      auto& first_child = localView.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      assert(((fe.size()== 4 and Traits::mydim==2) or (fe.size()== 8 and Traits::mydim==3)) && "EAS only supported for Q1 or H1 elements");
      for (int i = 0; i < fe.size(); ++i) {
        const int I    = Traits::worlddim * i;
        for (int j = 0; j < fe.size(); ++j) {
          const int J    = Traits::worlddim * j;
          Eigen::Matrix<double,enhancedStrainSize,strainSize> L;
          Eigen::Matrix<double,enhancedStrainSize,enhancedStrainSize> D;
          L.setZero();
          D.setZero();

          for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
            const auto Jinv = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().inverse().eval();
            const auto Bi = strainFunction.evaluateDerivative(gpIndex,wrt(coeff(J)), transformWith(Jinv));
            const auto Ceval = C(gpIndex);
            const auto M = easFunction.calcM(gp.position());
            const double detJ = geo.integrationElement(gp.position());

            L += M.transpose() * Ceval * Bi * detJ * gp.weight();
            D += M.transpose() * Ceval * M * detJ * gp.weight();
          }
          h.template block<Traits::worldDim,Traits::worldDim>(I,J)-= L.transpose()*D.inverse()*L;
        }
      }

    }

void setEASType(EASType otherEASType)
    {
  easType =otherEASType;
}

private:
  DisplacementBasedElement displacementBasedElement;
  EASType easType{EASType::Q1E4};
};




}
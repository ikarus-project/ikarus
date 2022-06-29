
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

namespace Ikarus{
class EnhancedAssumedStrains {


  virtual void transformStiffness(K,Bop, const Geometry& geo)=0;
  void transformInternalForces(K,Bop)=0;


};

class EAS_E4 : public EnhancedAssumedStrains
{

  void transformStiffness(K,Bop) final
  {
    Dune::FieldVector<double,2> quadPos0;
    quadPos0[0] = 0.5; // Center of the Element in Domain [0,1]
    quadPos0[1] = 0.5; // Center of the Element in Domain [0,1]

    const auto jacobianinvT0 = geometry.jacobianInverseTransposed(quadPos0); //J^{-1}.Transpose() in Dune = J^{-1}
    const auto detJ0 = geometry.integrationElement(quadPos0); //determinant(J)

    Eigen::MatrixXd jaco_it = Eigen::MatrixXd::Zero(2,2);
    for (size_t i=0;i<2;i++)
      for (size_t j=0;j<2;j++)
        jaco_it(i,j)=jacobianinvT0[i][j];

    auto jaco = (jaco_it).inverse();
    auto J11 = jaco(0,0);
    auto J12 = jaco(0,1);
    auto J21 = jaco(1,0);
    auto J22 = jaco(1,1);

    Eigen::Matrix3d T0;
    T0 << J11*J11 , J12*J12 , J11*J12 ,
        J21*J21 , J22*J22 , J21*J22 ,
        2.0*J11*J21 , 2.0*J12*J22 , J21*J12 + J11*J22;
    auto T0_inv = T0.inverse() * (detJ0/detJ);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(3, 4);
    M(0,0) = quadPos[0]-0.5;
    M(1,1) = quadPos[1]-0.5;
    M(2,2) = quadPos[0]-0.5;
    M(2,3) = quadPos[1]-0.5;

    M = T0_inv * M ;
  }
};


}
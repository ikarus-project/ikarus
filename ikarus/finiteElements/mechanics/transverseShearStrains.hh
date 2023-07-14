// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <Eigen/Core>
#include <dune/common/fvector.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/common/overloadset.hh>
#include <ranges>
namespace Ikarus {



struct DefaultTransverseShear {

  template<typename Geometry>
  void pre( const Geometry &geo,
            const auto &uFunction){}

  template< typename Geometry>
  auto value(const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex,
             const Geometry &geo,const auto& kin,
             const auto &uFunction, const auto& directorFunction, const auto& referenceDirectorFunction) const -> Eigen::Vector2<typename std::remove_cvref_t<decltype(uFunction)>::ctype>{
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    const auto J = Dune::toEigen(geo.jacobianTransposed(gpPos));
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
        uFunction.evaluateDerivative(gpPos, //Here the gpIndex could be passed
                                     Dune::wrt(spatialAll),
                                     Dune::on(Dune::DerivativeDirections::referenceElement)));
    const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();
    Eigen::Vector2<ScalarType> gammaV;
    const Eigen::Vector<ScalarType, 3> t0
        = referenceDirectorFunction.evaluate(gpPos,Dune::on(Dune::DerivativeDirections::referenceElement));
    const Eigen::Vector<ScalarType, 3> t
        = directorFunction.evaluate(gpPos,Dune::on(Dune::DerivativeDirections::referenceElement));
    gammaV<<j.row(0).dot(t) - J.row(0).dot(t0),j.row(1).dot(t) - J.row(1).dot(t0);

    return gammaV;
  }
  template<typename Geometry>
  auto derivativeWRTMidSurface(const auto& kin,const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex,
                               const auto &dNAtGp, const Geometry& geo,const auto& uFunction,const auto& directorFunction,
                               const auto& localBasis,                  const int node)  const{
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

    using namespace Dune::TypeTree::Indices;
    using namespace Dune::DerivativeDirections;
    const Eigen::Vector<ScalarType, 3> t
        = directorFunction.evaluate(gpPos,Dune::on(Dune::DerivativeDirections::referenceElement));
    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda2
        = uFunction.evaluateDerivative(gpPos, Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(Dune::DerivativeDirections::referenceElement));
    const auto &diffa1 = diffa1Anda2[0];
    const auto &diffa2 = diffa1Anda2[1];
    Eigen::Matrix<ScalarType, 2, 3> bop;
    bop.setZero();
    bop.row(0) = t.transpose() * diffa1(0,0);  // trans_shear_{,disp}
    bop.row(1) =t.transpose() * diffa2(0,0);
    return bop;
  }

  template<typename Geometry>
  auto derivativeWRTDirector(const auto& kin,const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex,
                             const auto &dNAtGp, const Geometry& geo,const auto& uFunction,const auto& directorFunction, const auto& localBasis,
                               const int node)  const{
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

    using namespace Dune::TypeTree::Indices;
    using namespace Dune::DerivativeDirections;
    const Eigen::Matrix<ScalarType, 3, 2> difft
        = directorFunction.evaluateDerivative(gpPos, Dune::wrt(coeff(_1, node)),Dune::on(Dune::DerivativeDirections::referenceElement));

    Eigen::Matrix<ScalarType, 2, 2> bop;
    const auto J = Dune::toEigen(geo.jacobianTransposed(gpPos));
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
        uFunction.evaluateDerivative(gpPos, //Here the gpIndex could be passed
                                     Dune::wrt(spatialAll),
                                     Dune::on(Dune::DerivativeDirections::referenceElement)));
    const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();

    bop.row(0) =j.row(0) * difft;  // trans_shear_{,disp}
    bop.row(1) = j.row(1) * difft;

    return bop;
  }


  template<typename Geometry,typename ScalarType>
  auto secondDerivativeWRTDirectorDirector(const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex,
                                           const auto &dNAtGp, const Geometry& geo,const auto& uFunction,
                                           const auto& directorFunction, const auto& localBasis,
                        const Eigen::Vector<ScalarType,8> &S,const auto& kin, int i, int j)const {
    using namespace Dune::TypeTree::Indices;
//    using namespace Dune::DerivativeDirections;
//    const Eigen::Matrix<ScalarType, 3, 2> difft
//        = directorFunction.evaluateDerivative(gpPos, Dune::wrt(coeff(_1, node)),Dune::on(Dune::DerivativeDirections::referenceElement));

    Eigen::Matrix<ScalarType, 2, 2> bop;
    const auto J = Dune::toEigen(geo.jacobianTransposed(gpPos));
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
        uFunction.evaluateDerivative(gpPos, //Here the gpIndex could be passed
                                     Dune::wrt(spatialAll),
                                     Dune::on(Dune::DerivativeDirections::referenceElement)));
    const Eigen::Matrix<ScalarType, 2, 3> jE = J + gradu.transpose();

    const auto a1 = jE.row(0).transpose().eval();
    const auto a2 = jE.row(1).transpose().eval();
    const Eigen::Matrix<ScalarType, 2, 2> S1
        = directorFunction.evaluateDerivative(gpPos, Dune::wrt(coeff(_1, i,_1, j)),Dune::along(a1),Dune::on(Dune::DerivativeDirections::referenceElement));

    const Eigen::Matrix<ScalarType, 2, 2> S2
        = directorFunction.evaluateDerivative(gpPos, Dune::wrt(coeff(_1, i,_1, j)),Dune::along(a2),Dune::on(Dune::DerivativeDirections::referenceElement));
    Eigen::Matrix<ScalarType, 2, 2> kg = S1* S[6] + S2 * S[7];

    return kg;
  }

  template<typename Geometry,typename ScalarType>
  auto secondDerivativeWRTSurfaceDirector(const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex,
                                          const auto &dNAtGp, const Geometry& geo,const auto& uFunction,
                                          const auto& directorFunction, const auto& localBasis,
                                           const Eigen::Vector<ScalarType,8> &S,const auto& kin, int i, int j)const {
    using namespace Dune::TypeTree::Indices;
    using namespace Dune::DerivativeDirections;
    const double& dN1i = dNAtGp(i, 0);
    const double& dN2i = dNAtGp(i, 1);
    const Eigen::Matrix<ScalarType, 3, 2> P
        = directorFunction.evaluateDerivative(gpPos, Dune::wrt(coeff(_1, j)),Dune::on(Dune::DerivativeDirections::referenceElement));
    Eigen::Matrix<ScalarType, 3, 2> kg = P * (dN1i * S[6] + dN2i * S[7]); // shear_{,dir,disp}*Q
    // bending_{,dir,disp}*M

    return kg;

    return kg;
  }

};


struct CASAnsatzFunctionTShear
{
  Dune::LagrangeCubeLocalFiniteElement<double, double, 2, 1> q1lfem2D;
  mutable std::vector<double> out;

  auto positions(std::vector<Dune::FieldVector<double, 2>>& lagrangePoints)
  {
    lagrangePoints.resize(q1lfem2D.size());
    for (int i = 0; i < 2; i++) {
      auto ithCoord = [&i](const Dune::FieldVector<double, 2>& x) { return x[i]; };
      q1lfem2D.localInterpolation().interpolate(ithCoord, out);
      for (std::size_t jI = 0; jI < out.size(); jI++)
        lagrangePoints[jI][i] = out[jI];
    }
  }

  void evaluateFunction(const Dune::FieldVector<double,2>& gpPos,std::vector<Dune::FieldVector<double, 1>>& N) const
  {
    q1lfem2D.localBasis().evaluateFunction(gpPos, N);
  }
};

struct CASAnsatzFunctionANSTShear
{
  Dune::LagrangeCubeLocalFiniteElement<double, double, 2, 1> q1lfem2D;
  mutable std::vector<double> out;

  auto positions(std::vector<Dune::FieldVector<double, 2>>& lagrangePoints)
  {
    lagrangePoints.resize(4);
    lagrangePoints[0] = {0.5,0};
    lagrangePoints[1] = {0,0.5};
    lagrangePoints[2] = {0.5, 1};
    lagrangePoints[3] = {1, 0.5};
  }

  void evaluateFunction(const Dune::FieldVector<double,2>& gpPos,std::vector<Dune::FieldVector<double, 1>>& NANS) const
  {
    NANS.resize(4);
    NANS[0] =  (1 - gpPos[1]);
    NANS[1] = (1 - gpPos[0]);
    NANS[2] = gpPos[1];
    NANS[3] = gpPos[0];
  }

};

struct CASTransverseStrain
{
  DefaultTransverseShear defaultTransverseShear;
  mutable std::vector<Dune::FieldVector<double, 2>> lagrangePoints;
  CASAnsatzFunctionANSTShear cASAnsatzFunction;


  template<typename T>
  using Vec = std::vector<Eigen::Vector<T, 3>>;
  std::variant<Vec<double>,Vec<autodiff::dual>,Vec<autodiff::dual2nd>> membraneStrainsAtVertices;
  mutable std::vector<Dune::FieldVector<double, 1>> NANS;

  CASTransverseStrain()
  {
    cASAnsatzFunction.positions(lagrangePoints);
}
  template<typename Geometry>
  void pre( const Geometry &geo,
           const auto &uFunction) {
//    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

//    using namespace Dune::DerivativeDirections;
//    using namespace Dune;
//
//    membraneStrainsAtVertices = std::vector<Eigen::Vector<ScalarType, 3>>();
//    for (int i = 0; auto& lP : lagrangePoints) {
//      const auto J                                = toEigen(geo.jacobianTransposed(lP));
//      const Eigen::Matrix<double, 2, 2> A         = J * J.transpose();
//      const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
//        uFunction.evaluateDerivative(lP, wrt(spatialAll, Dune::on(DerivativeDirections::referenceElement))));
//
//      const auto  epsV = defaultMembraneStrain.value(lP,geo,uFunction);
//
//      std::visit([&](auto& vec){
//        if constexpr (std::is_same_v<typename std::remove_cvref_t<decltype(vec)>::value_type,ScalarType>) {
//          vec.push_back(epsV);}
//          },membraneStrainsAtVertices);
//    }

  }
  template< typename Geometry>
  auto value(const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex,
             const Geometry &geo,const auto& kin,
             const auto &uFunction, const auto& directorFunction, const auto& referenceDirectorFunction) const -> Eigen::Vector2<typename std::remove_cvref_t<decltype(uFunction)>::ctype>{
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

    cASAnsatzFunction.evaluateFunction(gpPos, NANS);

    Eigen::Vector<ScalarType, 2> res;
    res.setZero();
    if(lagrangePoints.size()!=4)
      DUNE_THROW(Dune::Exception,"Wrong node size");
   res[0]= NANS[0][0]*defaultTransverseShear.value(lagrangePoints[0],integrationPointIndex,geo,kin,uFunction,directorFunction,referenceDirectorFunction)[0]+
    NANS[2][0]*defaultTransverseShear.value(lagrangePoints[2],integrationPointIndex,geo,kin,uFunction,directorFunction,referenceDirectorFunction)[0];

   res[1]= NANS[1][0]*defaultTransverseShear.value(lagrangePoints[1],integrationPointIndex,geo,kin,uFunction,directorFunction,referenceDirectorFunction)[1]+
            NANS[3][0]*defaultTransverseShear.value(lagrangePoints[3],integrationPointIndex,geo,kin,uFunction,directorFunction,referenceDirectorFunction)[1];
//    std::visit([&](auto& vec){
//      if constexpr (std::is_same_v<typename std::remove_cvref_t<decltype(vec)>::value_type,ScalarType>) {
//        for (int i = 0; i < NANS.size(); ++i) {
//          res += vec[i] * NANS[i][0];
//        }
//      }
//      ;},membraneStrainsAtVertices);

    return res;
  }
  mutable Eigen::Matrix<double,Eigen::Dynamic,2> dN;
  template<typename Geometry>
  auto derivativeWRTMidSurface(const auto& kin,const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex, const auto &dNAtGp, const Geometry& geo,const auto& uFunction,const auto& directorFunction, const auto& localBasis,
                               const int node)  const{
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

    cASAnsatzFunction.evaluateFunction(gpPos, NANS);
    Eigen::Matrix<ScalarType, 2, 3> bop;
    bop.setZero();
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    using namespace Dune::Indices;
//    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda20
//        = uFunction.evaluateDerivative(lagrangePoints[0], Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(referenceElement));
//    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda21
//        = uFunction.evaluateDerivative(lagrangePoints[1], Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(referenceElement));
//    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda21
//        = uFunction.evaluateDerivative(lagrangePoints[2], Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(referenceElement));
//    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda21
//        = uFunction.evaluateDerivative(lagrangePoints[3], Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(referenceElement));

    bop.row(0) = NANS[0][0]*defaultTransverseShear.derivativeWRTMidSurface(kin,lagrangePoints[0],integrationPointIndex,dNAtGp,geo,uFunction,directorFunction,localBasis,node).row(0)+
             NANS[2][0]*defaultTransverseShear.derivativeWRTMidSurface(kin,lagrangePoints[2],integrationPointIndex,dNAtGp,geo,uFunction,directorFunction,localBasis,node).row(0);

    bop.row(1) = NANS[1][0]*defaultTransverseShear.derivativeWRTMidSurface(kin,lagrangePoints[1],integrationPointIndex,dNAtGp,geo,uFunction,directorFunction,localBasis,node).row(1)+
                 NANS[3][0]*defaultTransverseShear.derivativeWRTMidSurface(kin,lagrangePoints[3],integrationPointIndex,dNAtGp,geo,uFunction,directorFunction,localBasis,node).row(1);
//    for (int i = 0; auto& lP : lagrangePoints) {
//       localBasis.evaluateJacobian(lP,dN);
////       std::cout<<"lP: "<<lP<<std::endl;
////       std::cout<<dN<<std::endl;
//
//
//      const auto J = toEigen(geo.jacobianTransposed(lP));
//      const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
//        uFunction.evaluateDerivative(lP, wrt(spatialAll), Dune::on(DerivativeDirections::referenceElement)));
//      const Eigen::Matrix<ScalarType, 2, 3> jE = J + gradu.transpose();
//      const auto  bopI = defaultTransverseShear.derivative<Geometry,ScalarType>(lP,jE,dN,geo,uFunction,localBasis,node);
//      bop+=bopI*NANS[i++][0];
//    }
    return bop;
  }

  template<typename Geometry>
  auto derivativeWRTDirector(const auto& kin,const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex, const auto &dNAtGp, const Geometry& geo,const auto& uFunction,const auto& directorFunction, const auto& localBasis,
                               const int node)  const{
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

    cASAnsatzFunction.evaluateFunction(gpPos, NANS);
    Eigen::Matrix<ScalarType, 2, 2> bop;
    bop.setZero();
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
//    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda20
//        = uFunction.evaluateDerivative(lagrangePoints[0], Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(referenceElement));
//    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda21
//        = uFunction.evaluateDerivative(lagrangePoints[1], Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(referenceElement));
//    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda21
//        = uFunction.evaluateDerivative(lagrangePoints[2], Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(referenceElement));
//    const std::array<Eigen::Matrix<ScalarType, 3, 3>, 2> diffa1Anda21
//        = uFunction.evaluateDerivative(lagrangePoints[3], Dune::wrt(spatialAll, coeff(_0, node)),Dune::on(referenceElement));

    bop.row(0) = NANS[0][0]*defaultTransverseShear.derivativeWRTDirector(kin,lagrangePoints[0],integrationPointIndex,dNAtGp,geo,uFunction,directorFunction,localBasis,node).row(0)+
                 NANS[2][0]*defaultTransverseShear.derivativeWRTDirector(kin,lagrangePoints[2],integrationPointIndex,dNAtGp,geo,uFunction,directorFunction,localBasis,node).row(0);

    bop.row(1) = NANS[1][0]*defaultTransverseShear.derivativeWRTDirector(kin,lagrangePoints[1],integrationPointIndex,dNAtGp,geo,uFunction,directorFunction,localBasis,node).row(1)+
                 NANS[3][0]*defaultTransverseShear.derivativeWRTDirector(kin,lagrangePoints[3],integrationPointIndex,dNAtGp,geo,uFunction,directorFunction,localBasis,node).row(1);
    //    for (int i = 0; auto& lP : lagrangePoints) {
    //       localBasis.evaluateJacobian(lP,dN);
    ////       std::cout<<"lP: "<<lP<<std::endl;
    ////       std::cout<<dN<<std::endl;
    //
    //
    //      const auto J = toEigen(geo.jacobianTransposed(lP));
    //      const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
    //        uFunction.evaluateDerivative(lP, wrt(spatialAll), Dune::on(DerivativeDirections::referenceElement)));
    //      const Eigen::Matrix<ScalarType, 2, 3> jE = J + gradu.transpose();
    //      const auto  bopI = defaultTransverseShear.derivative<Geometry,ScalarType>(lP,jE,dN,geo,uFunction,localBasis,node);
    //      bop+=bopI*NANS[i++][0];
    //    }
    return bop;
  }

  template<typename Geometry,typename ScalarType>
  auto secondDerivativeWRTDirectorDirector(const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex,
                                           const auto &dNAtGp, const Geometry& geo,const auto& uFunction,
                                           const auto& directorFunction, const auto& localBasis,
                                           const  Eigen::Vector<ScalarType,8> &S,const auto& kin, int i, int j)const {
    cASAnsatzFunction.evaluateFunction(gpPos, NANS);
    Eigen::Matrix<ScalarType, 2, 2> kg;
    kg.setZero();
    using namespace Dune::DerivativeDirections;
    using namespace Dune;

    Eigen::Vector<ScalarType,8> S0;
    S0<<0,0,0,0,0,0,S[6],0;
    localBasis.evaluateJacobian(lagrangePoints[0],dN);
    const auto  kgIJ0 = defaultTransverseShear.secondDerivativeWRTDirectorDirector(lagrangePoints[0],integrationPointIndex,dN,geo,uFunction,directorFunction,localBasis,S0,kin,i,j);
    localBasis.evaluateJacobian(lagrangePoints[2],dN);
    const auto  kgIJ2 = defaultTransverseShear.secondDerivativeWRTDirectorDirector(lagrangePoints[2],integrationPointIndex,dN,geo,uFunction,directorFunction,localBasis,S0,kin,i,j);
    kg=NANS[0][0]*kgIJ0+NANS[2][0]*kgIJ2;
    S0<<0,0,0,0,0,0,0,S[7];
    localBasis.evaluateJacobian(lagrangePoints[1],dN);
    const auto  kgIJ1 = defaultTransverseShear.secondDerivativeWRTDirectorDirector(lagrangePoints[1],integrationPointIndex,dN,geo,uFunction,directorFunction,localBasis,S0,kin,i,j);
    localBasis.evaluateJacobian(lagrangePoints[3],dN);
    const auto  kgIJ3 = defaultTransverseShear.secondDerivativeWRTDirectorDirector(lagrangePoints[3],integrationPointIndex,dN,geo,uFunction,directorFunction,localBasis,S0,kin,i,j);
    kg+=NANS[1][0]*kgIJ1+NANS[3][0]*kgIJ3;

    return kg;
  }

  template<typename Geometry,typename ScalarType>
  auto secondDerivativeWRTSurfaceDirector(const Dune::FieldVector<double, 2> &gpPos,const int integrationPointIndex,
                                           const auto &dNAtGp, const Geometry& geo,const auto& uFunction,
                                           const auto& directorFunction, const auto& localBasis,
                                           const Eigen::Vector<ScalarType,8> &S,const auto& kin, int i, int j)const {
    cASAnsatzFunction.evaluateFunction(gpPos, NANS);
    Eigen::Matrix<ScalarType, 3, 2> kg;
    kg.setZero();
    using namespace Dune::DerivativeDirections;
    using namespace Dune;

    Eigen::Vector<ScalarType,8> S0;
    S0<<0,0,0,0,0,0,S[6],0;
    localBasis.evaluateJacobian(lagrangePoints[0],dN);
    const auto  kgIJ0 = defaultTransverseShear.secondDerivativeWRTSurfaceDirector(lagrangePoints[0],integrationPointIndex,dN,geo,uFunction,directorFunction,localBasis,S0,kin,i,j);
    localBasis.evaluateJacobian(lagrangePoints[2],dN);
    const auto  kgIJ2 = defaultTransverseShear.secondDerivativeWRTSurfaceDirector(lagrangePoints[2],integrationPointIndex,dN,geo,uFunction,directorFunction,localBasis,S0,kin,i,j);
    kg=NANS[0][0]*kgIJ0+NANS[2][0]*kgIJ2;
    S0<<0,0,0,0,0,0,0,S[7];
    localBasis.evaluateJacobian(lagrangePoints[1],dN);
    const auto  kgIJ1 = defaultTransverseShear.secondDerivativeWRTSurfaceDirector(lagrangePoints[1],integrationPointIndex,dN,geo,uFunction,directorFunction,localBasis,S0,kin,i,j);
    localBasis.evaluateJacobian(lagrangePoints[3],dN);
    const auto  kgIJ3 = defaultTransverseShear.secondDerivativeWRTSurfaceDirector(lagrangePoints[3],integrationPointIndex,dN,geo,uFunction,directorFunction,localBasis,S0,kin,i,j);
    kg+=NANS[1][0]*kgIJ1+NANS[3][0]*kgIJ3;
    return kg;
  }

};


template<class... Implementations>
class TransverseShearStrainVariant
{

 public:

  template<class Implementation>
  explicit TransverseShearStrainVariant(const Implementation& impl) :
      impl_(impl)
  {}

  TransverseShearStrainVariant() = default;
  TransverseShearStrainVariant(const TransverseShearStrainVariant& other) = default;
  template<class Implementation> requires (!std::is_same_v<Implementation,TransverseShearStrainVariant>)
  TransverseShearStrainVariant& operator=(const Implementation& impl)
      {
        impl_=impl;
        return *this;
      };
  TransverseShearStrainVariant(TransverseShearStrainVariant&& other)  noexcept = default;
  TransverseShearStrainVariant& operator=(const TransverseShearStrainVariant& other) = default;
  TransverseShearStrainVariant& operator=(TransverseShearStrainVariant&& other)  noexcept = default;

  template<typename Geometry>
  void pre( const Geometry &geo,
            const auto &uFunction){
  std::visit([&]( auto& impl) { impl.pre(geo,uFunction); }, impl_);
  }
  template<typename... Args>
  auto value(Args&&... args)  {

  return std::visit([&](const auto& impl) { return impl.value(std::forward<Args>(args)...); }, impl_);
  }
  template<typename... Args>
  auto secondDerivativeWRTDirectorDirector(Args&&... args)  {

    return std::visit([&](const auto& impl) { return impl.secondDerivativeWRTDirectorDirector(std::forward<Args>(args)...); }, impl_);
  }

  template<typename... Args>
  auto secondDerivativeWRTSurfaceDirector(Args&&... args)  {

    return std::visit([&](const auto& impl) { return impl.secondDerivativeWRTSurfaceDirector(std::forward<Args>(args)...); }, impl_);
  }

  template<typename... Args>
  auto derivativeWRTMidSurface(Args&&... args)  {

    return std::visit([&](const auto& impl) { return impl.derivativeWRTMidSurface(std::forward<Args>(args)...); }, impl_);
  }
  template<typename... Args>
  auto derivativeWRTDirector(Args&&... args)  {

    return std::visit([&](const auto& impl) { return impl.derivativeWRTDirector(std::forward<Args>(args)...); }, impl_);
  }

 private:
  std::variant<Implementations...> impl_;
};

using TransverseShearStrain =TransverseShearStrainVariant<DefaultTransverseShear,CASTransverseStrain>;

}

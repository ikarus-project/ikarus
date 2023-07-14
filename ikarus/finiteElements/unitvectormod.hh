//
// Created by Alex on 19.05.2021.
//

#pragma once
#include <concepts>

#include <ikarus/utils/eigenDuneTransformations.hh>
#include <dune/localfefunctions/manifolds/unitVector.hh>

namespace Dune {
const double unitvectol = 1e-8;
  template<typename Scalar>
  Scalar arcCosSquared(const Scalar& x) {
    using std::acos;
    const Scalar eps = 1e-2;
    if (x > 1-eps) {  // acos is not differentiable, use the series expansion instead,
      // we need here lots of terms to be sure that the numerical derivatives are also within maschine precission
      return 11665028.0/4729725.0
             -141088.0/45045.0*x
             +   413.0/429.0*x*x
             -  5344.0/12285.0*Dune::power(x,3)
             +    245.0/1287.0*Dune::power(x,4)
             -  1632.0/25025.0*Dune::power(x,5)
             +     56.0/3861.0*Dune::power(x,6)
             -    32.0/21021.0*Dune::power(x,7);
    } else {
      return Dune::power(acos(x),2);
    }
  }

  template<typename Scalar>
  Scalar derivativeOfArcCosSquared(const Scalar& x) {
    using std::acos;
    using std::sqrt;
    const Scalar eps = 1e-2;
    if (x > 1-eps) {  // regular expression is unstable, use the series expansion instead
      // we need here lots of terms to be sure that the numerical derivatives are also within maschine precission
      //return -2 + 2*(x-1)/3 - 4/15*(x-1)*(x-1);
      return -47104.0/15015.0
             +12614.0/6435.0*x
             -63488.0/45045.0*x*x
             + 1204.0/1287.0*Dune::power(x,3)
             - 2048.0/4095.0*Dune::power(x,4)
             +   112.0/585.0*Dune::power(x,5)
             -2048.0/45045.0*Dune::power(x,6)
             +   32.0/6435.0*Dune::power(x,7);

    } else if (x < -1+unitvectol) {  // The function is not differentiable
      DUNE_THROW(Dune::Exception, "arccos^2 is not differentiable at x==-1!"<< x);
    } else
      return -2*acos(x) / sqrt(1-x*x);
  }


  /** \brief Compute the second derivative of arccos^2 without getting unstable for x close to 1 */
  template<typename Scalar>
  Scalar secondDerivativeOfArcCosSquared(const Scalar& x) {
    using std::acos;
    using std::pow;
    const Scalar eps = 1e-2;
    if (x > 1-eps) {  // regular expression is unstable, use the series expansion instead
      // we need here lots of terms to be sure that the numerical derivatives are also within maschine precission
      //return 2.0/3 - 8*(x-1)/15;
      return 1350030.0/676039.0+5632.0/2028117.0*Dune::power(x,10)
             -1039056896.0/334639305.0*x
             +150876.0/39767.0*x*x
             -445186048.0/111546435.0*Dune::power(x,3)
             +       343728.0/96577.0*Dune::power(x,4)
             -  57769984.0/22309287.0*Dune::power(x,5)
             +      710688.0/482885.0*Dune::power(x,6)
             -  41615360.0/66927861.0*Dune::power(x,7)
             +     616704.0/3380195.0*Dune::power(x,8)
             -     245760.0/7436429.0*Dune::power(x,9);
    } else if (x < -1+unitvectol) {  // The function is not differentiable
      DUNE_THROW(Dune::Exception, "arccos^2 is not differentiable at x==-1!");
    } else
      return 2/(1-x*x) - 2*x*acos(x) / pow(1-x*x,1.5);
  }

  template<typename Scalar, int d>
  auto euclideanDerivativeOfDistanceSquaredWRTSecondArgument(const UnitVector<Scalar,d> & l, const UnitVector<Scalar,d> & r)  {
    const auto dotprod = l.getValue().dot(r.getValue());

    return (l.getValue()*derivativeOfArcCosSquared(dotprod)).eval();
  }

  template<typename Scalar, int d>
  auto derivativeOfDistanceSquaredWRTSecondArgument(const UnitVector<Scalar,d> & l, const UnitVector<Scalar,d> & r)  {
    return (l.orthonormalFrame().transpose()*euclideanDerivativeOfDistanceSquaredWRTSecondArgument(l,r)).eval();
  }

  template<typename Scalar, int d>
  Eigen::Matrix<Scalar,d-1,d-1> secondDerivativeOfDistanceSquaredWRTSecondArgument(const UnitVector<Scalar,d> & l, const UnitVector<Scalar,d> & r)  {

    const auto dotprod = l.getValue().dot(r.getValue());
    const auto BLA = r.orthonormalFrame();
    const auto gradEuk = euclideanDerivativeOfDistanceSquaredWRTSecondArgument(l,r);
    return (BLA.transpose()*secondDerivativeOfArcCosSquared(dotprod)*l.getValue()*l.getValue().transpose()*BLA-r.getValue().dot(gradEuk)*Eigen::Matrix<double,d-1,d-1>::Identity()).eval();
  }

  template<typename Scalar, int d>
  Eigen::Matrix<Scalar,d,d> euclideanSecondDerivativeOfDistanceSquaredWRTSecondArgument(const UnitVector<Scalar,d> & l, const UnitVector<Scalar,d> & r)  {

    const auto dotprod = l.getValue().dot(r.getValue());
    const auto gradEuk = euclideanDerivativeOfDistanceSquaredWRTSecondArgument(l,r);
    return (secondDerivativeOfArcCosSquared(dotprod)*l.getValue()*l.getValue().transpose()-r.getValue().dot(gradEuk)*Eigen::Matrix<double,d,d>::Identity()).eval();
  }
}  // namespace Ikarus

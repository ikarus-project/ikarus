//
// Created by Alex on 12.05.2021.
//



#ifndef IKARUS_INTERPOLATOR_H
#define IKARUS_INTERPOLATOR_H

#include <numeric>

template<typename DerivedAnsatzFunctionType,typename DerivedCoeffsType>
auto  interpolate(const Eigen::MatrixBase<DerivedAnsatzFunctionType>& N ,const Eigen::MatrixBase<DerivedCoeffsType>& coeffs)
{
        return (coeffs*N).eval();
}

//template<typename DerivedAnsatzFunctionType,typename GlobalCoordinateListType>
//auto  interpolateRowWise(const std::span<AnsatzfunctionType> ansatzfunction ,const std::span<CoeffsType>& coeffs)
//{
//    assert(ansatzfunction.size()==coeffs.size());
//    if constexpr (std::is_same_v<double,AnsatzfunctionType>)
//        return std::inner_product(ansatzfunction.begin(),ansatzfunction.end(),coeffs.begin(),typename CoeffsType::CoordinateType {},
//                                  std::plus<typename CoeffsType::CoordinateType>(),[&](const double& a, const CoeffsType& b){
//                    return a*(b);
//                });
//
//}

#endif //IKARUS_INTERPOLATOR_H

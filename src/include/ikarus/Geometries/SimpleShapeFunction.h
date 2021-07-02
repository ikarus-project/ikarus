//
// Created by ac120950 on 17.07.2020.
//

#pragma once
//#include "../Interfaces/InterfaceShapeFunction.h"
#include <functional>

// class SimpleShapeFunction {
//
// public:
//    typedef Eigen::Vector2d CoordinateType;
//    typedef double ctype;
//     typedef std::array<const std::function<double(const CoordinateType&) >,4> FunctionVectorType;
//     typedef Eigen::Matrix<ctype,4,1> VectorType;
//     typedef Eigen::Matrix<ctype,2,4> MatrixType;
//
//
//    SimpleShapeFunction(const FunctionVectorType & myshapeFunctions) :
//    shapeFunctions{myshapeFunctions}
//    {}
//
//    VectorType evaluate(CoordinateType coord) const
//    {
//        VectorType res;
//        res.resize(shapeFunctions.size());
//        int counter=0;
//        for (const auto& Ni:shapeFunctions) {
//            res(counter) = Ni(coord);
//            counter++;
//        }
//        return res;
//    }
//    MatrixType derivative(CoordinateType paraPoint) const{return MatrixType();};
//    int size() const{return shapeFunctions.size();};
//
// private:
//    FunctionVectorType shapeFunctions;
//};

//#endif  // FE_DESIGNPATTERN_SIMPLESHAPEFUNCTION_H

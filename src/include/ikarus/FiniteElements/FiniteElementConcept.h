//
// Created by Alex on 10.05.2021.
//

#pragma once
//#include <memory>
//
//#include "ikarus/utils/LinearAlgebraTypedefs.h"

// namespace Ikarus::Concepts {
//    template <typename FEConceptType>
//    concept FE = requires (FEConceptType fe){
//        typename FEConceptType::ctype;
//        FEConceptType::coorddimension;
//        FEConceptType::mydimension;
//        typename FEConceptType::Geometry;
//        typename FEConceptType::VectorType;
//        typename FEConceptType::MatrixType;
//        typename FEConceptType::DofVectorType;
//        { fe.initialize() } ->  std::same_as<void>;
//        { fe.dofSize() } ->  std::same_as<int>;
//        { fe.getEnergy() } ->  std::same_as<typename FEConceptType::ctype>;
//        { fe.calculateLocalSystem() } ->  std::same_as<std::pair<typename
//        FEConceptType::VectorType,
//                                                                                      typename
//                                                                                      FEConceptType::MatrixType>>;
//        { fe.calculateLHS() } ->  std::same_as<typename FEConceptType::MatrixType>;
//        { fe.calculateRHS() } ->  std::same_as<typename FEConceptType::VectorType>;
//        { fe.getDofVector() } ->  std::same_as<typename FEConceptType::DofVectorType>;
//    };

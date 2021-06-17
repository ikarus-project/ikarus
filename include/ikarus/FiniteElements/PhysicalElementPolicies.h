//
// Created by Alex on 16.06.2021.
//

#pragma once


namespace Ikarus::Impl::Concepts
{

#define TRYCALLMEMBERFUNCTION(Str) \
   if constexpr(Ikarus::Impl::Concepts::Has##Str<FE>) \
return fe.Str(); \
else \
DUNE_THROW(Dune::InvalidStateException, "The member function \"" << #Str << "\" is not implemented by this element"); \


#define TRYCALLMEMBERFUNCTIONDONTTHROW(Str) \
   if constexpr(Ikarus::Impl::Concepts::Has##Str<FE>) \
return fe.Str(); \
else return; \


template<typename PhysicalType>
concept HascalculateLHS = requires(PhysicalType pfe){
  pfe.calculateLHS();
};

template<typename PhysicalType>
concept HascalculateRHS = requires(PhysicalType pfe){
  pfe.calculateRHS();
};

template<typename PhysicalType>
concept HasgetDofVector = requires(PhysicalType pfe){
  pfe.getDofVector();
};


template<typename PhysicalType>
concept HasdofSize = requires(PhysicalType pfe){
  pfe.dofSize();
};

template<typename PhysicalType>
concept HascalculateLocalSystem = requires(PhysicalType pfe){
  pfe.calculateLocalSystem();
};

template<typename PhysicalType>
concept Hasinitialize = requires(PhysicalType pfe){
  pfe.initialize();
};
}
//
// Created by Alex on 16.06.2021.
//

#pragma once

namespace Ikarus::Concepts {

#define TRYCALLFUNCTION(Str)                             \
  if constexpr (Ikarus::Concepts::Has##Str<FE>)          \
    return fe.Str();                                     \
  else if constexpr (Ikarus::Concepts::HasFree##Str<FE>) \
    return Str(fe);                                      \
  else                                                   \
    DUNE_THROW(Dune::InvalidStateException,              \
               "The member function \"" << #Str << "\" is not implemented by this element");

#define TRYCALLFUNCTIONDONTTHROW(Str)                    \
  if constexpr (Ikarus::Concepts::Has##Str<FE>)          \
    return fe.Str();                                     \
  else if constexpr (Ikarus::Concepts::HasFree##Str<FE>) \
    return Str(fe);                                      \
  else                                                   \
    return;

  template <typename PhysicalType>
  concept HascalculateLHS = requires(PhysicalType pfe) {
    pfe.calculateLHS();
  };

  template <typename PhysicalType>
  concept HasFreecalculateLHS = requires(PhysicalType pfe) {
    calculateLHS(pfe);
  };

  template <typename PhysicalType>
  concept HascalculateRHS = requires(PhysicalType pfe) {
    pfe.calculateRHS();
  };

  template <typename PhysicalType>
  concept HasFreecalculateRHS = requires(PhysicalType pfe) {
    calculateRHS(pfe);
  };

  template <typename PhysicalType>
  concept HasgetDofVector = requires(PhysicalType pfe) {
    pfe.getDofVector();
  };

  template <typename PhysicalType>
  concept HasFreegetDofVector = requires(PhysicalType pfe) {
    getDofVector(pfe);
  };

  template <typename PhysicalType>
  concept HasdofSize = requires(PhysicalType pfe) {
    pfe.dofSize();
  };

  template <typename PhysicalType>
  concept HasFreedofSize = requires(PhysicalType pfe) {
    dofSize(pfe);
  };

  template <typename PhysicalType>
  concept HascalculateLocalSystem = requires(PhysicalType pfe) {
    pfe.calculateLocalSystem();
  };

  template <typename PhysicalType>
  concept HasFreecalculateLocalSystem = requires(PhysicalType pfe) {
    calculateLocalSystem(pfe);
  };

  template <typename PhysicalType>
  concept Hasinitialize = requires(PhysicalType pfe) {
    pfe.initialize();
  };

  template <typename PhysicalType>
  concept HasFreeinitialize = requires(PhysicalType pfe) {
    initialize(pfe);
  };
}  // namespace Ikarus::Concepts
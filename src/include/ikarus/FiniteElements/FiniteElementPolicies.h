//
// Created by Alex on 16.06.2021.
//

#pragma once

namespace Ikarus::Concepts {




#define TRYCALLFUNCTION(Str,...)                             \
  if constexpr (Ikarus::Concepts::Has##Str<FE>)          \
    return fe.Str(__VA_ARGS__);                                     \
  else if constexpr (Ikarus::Concepts::HasFree##Str<FE>) \
    return Str(fe, ## __VA_ARGS__);                                      \
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
  concept HascalculateMatrix = requires(PhysicalType pfe) {
    pfe.calculateLHS();
  };

  template <typename PhysicalType>
  concept HasFreecalculateMatrix = requires(PhysicalType pfe) {
    calculateLHS(pfe);
  };

  template <typename PhysicalType>
  concept HascalculateScalar = requires(PhysicalType pfe) {
    pfe.calculateScalar();
  };

  template <typename PhysicalType>
  concept HasFreecalculateScalar = requires(PhysicalType pfe) {
    calculateScalar(pfe);
  };

  template <typename PhysicalType>
  concept HascalculateVector = requires(PhysicalType pfe) {
    pfe.calculateRHS();
  };

  template <typename PhysicalType>
  concept HasFreecalculateVector = requires(PhysicalType pfe) {
    calculateRHS(pfe);
  };

  template <typename PhysicalType>
  concept HasgetEntityVariablePairs = requires(PhysicalType pfe) {
    pfe.getEntityVariablePairs();
  };

  template <typename PhysicalType>
  concept HasFreegetEntityVariablePairs = requires(PhysicalType pfe) {
    getEntityVariablePairs(pfe);
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


namespace Ikarus::FiniteElements
{
  enum class ElementVectorAffordances
  {
    internalforces
  };

  enum class ElementMatrixAffordances
  {
    stiffnessMatrix,
    dStiffnessMatrixdBucklingVector,
    massMatrix
  };

  enum class ElementScalarAffordances
  {
    potentialEnergy
  };

  inline constexpr ElementVectorAffordances internalforces = ElementVectorAffordances::internalforces;

  inline constexpr ElementMatrixAffordances stiffnessMatrix = ElementMatrixAffordances::stiffnessMatrix;
  inline constexpr ElementMatrixAffordances dStiffnessMatrixdBucklingVector = ElementMatrixAffordances::dStiffnessMatrixdBucklingVector;
  inline constexpr ElementMatrixAffordances massMatrix = ElementMatrixAffordances::massMatrix;

  inline constexpr ElementScalarAffordances potentialEnergy = ElementScalarAffordances::potentialEnergy;

}
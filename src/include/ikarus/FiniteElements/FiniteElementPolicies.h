//
// Created by Alex on 16.06.2021.
//

#pragma once

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#include <ikarus/FiniteElements/FiniteElementPolicies.h>
#endif

namespace Ikarus::FiniteElements {
  enum class ElementVectorAffordances { forces };

  enum class ElementMatrixAffordances { stiffness, stiffnessdiffBucklingVector, mass };

  enum class ElementScalarAffordances { potentialEnergy };

  inline constexpr ElementVectorAffordances forces = ElementVectorAffordances::forces;

  inline constexpr ElementMatrixAffordances stiffness = ElementMatrixAffordances::stiffness;
  inline constexpr ElementMatrixAffordances stiffnessdiffBucklingVector
      = ElementMatrixAffordances::stiffnessdiffBucklingVector;
  inline constexpr ElementMatrixAffordances mass = ElementMatrixAffordances::mass;

  inline constexpr ElementScalarAffordances potentialEnergy = ElementScalarAffordances::potentialEnergy;

}  // namespace Ikarus::FiniteElements

namespace Ikarus::Concepts {

#define TRYCALLFUNCTION(Str, ...)                        \
  if constexpr (Ikarus::Concepts::Has##Str<FE>)          \
    return fe.Str(__VA_ARGS__);                          \
  else if constexpr (Ikarus::Concepts::HasFree##Str<FE>) \
    return Str(fe, ##__VA_ARGS__);                       \
  else                                                   \
    DUNE_THROW(Dune::InvalidStateException,              \
               "The member/free function \"" << #Str << "\" is not implemented by this element");

#define TRYCALLFUNCTIONDONTTHROW(Str)                    \
  if constexpr (Ikarus::Concepts::Has##Str<FE>)          \
    return fe.Str();                                     \
  else if constexpr (Ikarus::Concepts::HasFree##Str<FE>) \
    return Str(fe);                                      \
  else                                                   \
    return;

  template <typename FiniteElement>
  concept HascalculateMatrix = requires(FiniteElement fe, Ikarus::FiniteElements::ElementMatrixAffordances matA) {
    fe.calculateMatrix(matA);
  };

  template <typename FiniteElement>
  concept HasFreecalculateMatrix = requires(FiniteElement fe, Ikarus::FiniteElements::ElementMatrixAffordances matA) {
    calculateMatrix(fe, matA);
  };

  template <typename FiniteElement>
  concept HascalculateScalar = requires(FiniteElement fe, Ikarus::FiniteElements::ElementScalarAffordances scalA) {
    fe.calculateScalar(scalA);
  };

  template <typename FiniteElement>
  concept HasFreecalculateScalar = requires(FiniteElement fe, Ikarus::FiniteElements::ElementScalarAffordances scalA) {
    calculateScalar(fe, scalA);
  };

  template <typename FiniteElement>
  concept HascalculateVector = requires(FiniteElement fe, Ikarus::FiniteElements::ElementVectorAffordances vecA) {
    fe.calculateVector(vecA);
  };

  template <typename FiniteElement>
  concept HasFreecalculateVector = requires(FiniteElement fe, Ikarus::FiniteElements::ElementVectorAffordances vecA) {
    calculateVector(fe, vecA);
  };

  template <typename FiniteElement>
  concept HasgetEntityVariablePairs = requires(FiniteElement fe) {
    fe.getEntityVariablePairs();
  };

  template <typename FiniteElement>
  concept HasFreegetEntityVariablePairs = requires(FiniteElement fe) {
    getEntityVariablePairs(fe);
  };

  template <typename FiniteElement>
  concept HasSomegetEntityVariablePairs
      = HasgetEntityVariablePairs<FiniteElement> || HasFreegetEntityVariablePairs<FiniteElement>;

  template <typename FiniteElement>
  concept HasdofSize = requires(FiniteElement fe) {
    fe.dofSize();
  };

  template <typename FiniteElement>
  concept HasFreedofSize = requires(FiniteElement fe) {
    dofSize(fe);
  };

  template <typename FiniteElement>
  concept HasgetEntityID = requires(FiniteElement fe) {
    fe.getEntityID();
  };

  template <typename FiniteElement>
  concept HasFreegetEntityID = requires(FiniteElement fe) {
    getEntityID(fe);
  };

  template <typename FiniteElement>
  concept HascalculateLocalSystem = requires(FiniteElement fe, Ikarus::FiniteElements::ElementMatrixAffordances matA,
                                             Ikarus::FiniteElements::ElementVectorAffordances vecA) {
    fe.calculateLocalSystem(matA, vecA);
  };

  template <typename FiniteElement>
  concept HasFreecalculateLocalSystem
      = requires(FiniteElement fe, Ikarus::FiniteElements::ElementMatrixAffordances matA,
                 Ikarus::FiniteElements::ElementVectorAffordances vecA) {
    calculateLocalSystem(fe, matA, vecA);
  };

  template <typename FiniteElement>
  concept Hasinitialize = requires(FiniteElement fe) {
    fe.initialize();
  };

  template <typename FiniteElement>
  concept HasFreeinitialize = requires(FiniteElement fe) {
    initialize(fe);
  };

  template <typename FiniteElement>
  concept MinimalFiniteElementLinearAlgebraAffordances
      = HasFreecalculateScalar<FiniteElement> || HasFreecalculateVector<FiniteElement> || HasFreecalculateMatrix<
          FiniteElement> || HascalculateScalar<FiniteElement> || HascalculateVector<FiniteElement> || HascalculateMatrix<FiniteElement>;

}  // namespace Ikarus::Concepts

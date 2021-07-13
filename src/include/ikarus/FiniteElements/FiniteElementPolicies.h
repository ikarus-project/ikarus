//
// Created by Alex on 16.06.2021.
//

#pragma once

//#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Variables/InterfaceVariable.h>
//#ifdef __clang__
//#  pragma clang diagnostic push
//#  pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
//#endif

namespace Ikarus::FiniteElements {
  enum class VectorAffordances { forces };

  enum class MatrixAffordances { stiffness, materialstiffness, geometricstiffness, stiffnessdiffBucklingVector, mass };

  enum class ScalarAffordances { potentialEnergy };

  inline constexpr VectorAffordances forces = VectorAffordances::forces;

  inline constexpr MatrixAffordances stiffness                   = MatrixAffordances::stiffness;
  inline constexpr MatrixAffordances stiffnessdiffBucklingVector = MatrixAffordances::stiffnessdiffBucklingVector;
  inline constexpr MatrixAffordances mass                        = MatrixAffordances::mass;

  inline constexpr ScalarAffordances potentialEnergy = ScalarAffordances::potentialEnergy;

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
  concept HascalculateMatrix = requires(FiniteElement fe, std::vector<Ikarus::Variable::IVariable*> vars,
                                        Ikarus::FiniteElements::MatrixAffordances matA) {
    fe.calculateMatrix(vars, matA);
  };

  template <typename FiniteElement>
  concept HasFreecalculateMatrix = requires(FiniteElement fe, std::vector<Ikarus::Variable::IVariable*> vars,
                                            Ikarus::FiniteElements::MatrixAffordances matA) {
    calculateMatrix(fe, vars, matA);
  };

  template <typename FiniteElement>
  concept HascalculateScalar = requires(FiniteElement fe, std::vector<Ikarus::Variable::IVariable*> vars,
                                        Ikarus::FiniteElements::ScalarAffordances scalA) {
    fe.calculateScalar(vars, scalA);
  };

  template <typename FiniteElement>
  concept HasFreecalculateScalar = requires(FiniteElement fe, std::vector<Ikarus::Variable::IVariable*> vars,
                                            Ikarus::FiniteElements::ScalarAffordances scalA) {
    calculateScalar(fe, vars, scalA);
  };

  template <typename FiniteElement>
  concept HascalculateVector = requires(FiniteElement fe, std::vector<Ikarus::Variable::IVariable*> vars,
                                        Ikarus::FiniteElements::VectorAffordances vecA) {
    fe.calculateVector(vars, vecA);
  };

  template <typename FiniteElement>
  concept HasFreecalculateVector = requires(FiniteElement fe, std::vector<Ikarus::Variable::IVariable*> vars,
                                            Ikarus::FiniteElements::VectorAffordances vecA) {
    calculateVector(fe, vars, vecA);
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
  concept HascalculateLocalSystem
      = requires(FiniteElement fe, std::vector<Ikarus::Variable::IVariable*> vars,
                 Ikarus::FiniteElements::MatrixAffordances matA, Ikarus::FiniteElements::VectorAffordances vecA) {
    fe.calculateLocalSystem(vars, matA, vecA);
  };

  template <typename FiniteElement>
  concept HasFreecalculateLocalSystem
      = requires(FiniteElement fe, std::vector<Ikarus::Variable::IVariable*> vars,
                 Ikarus::FiniteElements::MatrixAffordances matA, Ikarus::FiniteElements::VectorAffordances vecA) {
    calculateLocalSystem(fe, vars, matA, vecA);
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

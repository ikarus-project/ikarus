//
// Created by Alex on 16.06.2021.
//

#pragma once

//#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <optional>

#include <ikarus/Variables/InterfaceVariable.h>
#include <ikarus/Variables/ParameterFactory.h>
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#endif

namespace Ikarus {

  struct VariableIndicesPair {
    using VariableVector = std::vector<Variable::VariableTags>;
    using Indices        = std::vector<size_t>;
    Indices indices;
    VariableVector variableVector;
  };

  enum class VectorAffordances { noAffordance, forces };

  enum class MatrixAffordances {
    noAffordance,
    stiffness,
    materialstiffness,
    geometricstiffness,
    stiffnessdiffBucklingVector,
    mass
  };

  enum class ScalarAffordances { noAffordance, potentialEnergy };

  inline constexpr VectorAffordances forces = VectorAffordances::forces;

  inline constexpr MatrixAffordances stiffness                   = MatrixAffordances::stiffness;
  inline constexpr MatrixAffordances stiffnessdiffBucklingVector = MatrixAffordances::stiffnessdiffBucklingVector;
  inline constexpr MatrixAffordances mass                        = MatrixAffordances::mass;

  inline constexpr ScalarAffordances potentialEnergy = ScalarAffordances::potentialEnergy;

  template <typename SolutionVectorType=Eigen::VectorXd, typename ParameterType =std::map<decltype(FEParameterValuePair::type), double>>
  struct FErequirements {
    std::vector<std::reference_wrapper<const SolutionVectorType>> sols;
    ParameterType parameter;
    ScalarAffordances scalarAffordances{ScalarAffordances::noAffordance};
    VectorAffordances vectorAffordances{VectorAffordances::noAffordance};
    VectorAffordances matrixAffordances{MatrixAffordances::noAffordance};
  };

}  // namespace Ikarus::FiniteElements

namespace Ikarus::Concepts {

#define TRYCALLFUNCTIONANDRETURN(Str, ...)                         \
  if constexpr (Ikarus::Concepts::Has##Str<decltype(fe)>)          \
    return fe.Str(__VA_ARGS__);                                    \
  else if constexpr (Ikarus::Concepts::HasFree##Str<decltype(fe)>) \
    return Str(fe, ##__VA_ARGS__);                                 \
  else                                                             \
    DUNE_THROW(Dune::InvalidStateException,                        \
               "The member/free function \"" << #Str << "\" is not implemented by this element");

#define TRYCALLFUNCTIONDONTTHROW(Str)                              \
  if constexpr (Ikarus::Concepts::Has##Str<decltype(fe)>)          \
    return fe.Str();                                               \
  else if constexpr (Ikarus::Concepts::HasFree##Str<decltype(fe)>) \
    return Str(fe);                                                \
  else                                                             \
    return;

  template <typename FiniteElement>
  concept HascalculateMatrix = requires(FiniteElement fe, typename FiniteElement::FERequirementType req) {
    fe.calculateMatrix(req);
  };

  template <typename FiniteElement>
  concept HasFreecalculateMatrix = requires(FiniteElement fe, typename FiniteElement::FERequirementType req) {
    calculateMatrix(fe, req);
  };

  template <typename FiniteElement>
  concept HascalculateScalar = requires(FiniteElement fe, typename FiniteElement::FERequirementType req) {
    fe.calculateScalar(req);
  };

  template <typename FiniteElement>
  concept HasFreecalculateScalar = requires(FiniteElement fe, typename FiniteElement::FERequirementType req) {
    calculateScalar(fe, req);
  };

  template <typename FiniteElement>
  concept HascalculateVector = requires(FiniteElement fe, typename FiniteElement::FERequirementType req) {
    fe.calculateVector(req);
  };

  template <typename FiniteElement>
  concept HasFreecalculateVector = requires(FiniteElement fe, typename FiniteElement::FERequirementType req) {
    calculateVector(fe, req);
  };

  template <typename FiniteElement>
  concept HasglobalIndices = requires(FiniteElement fe) {
    fe.globalIndices();
  };

  template <typename FiniteElement>
  concept HasFreeglobalIndices = requires(FiniteElement fe) {
    globalIndices(fe);
  };

  template <typename FiniteElement>
  concept HasSomeglobalIndices = HasglobalIndices<FiniteElement> || HasFreeglobalIndices<FiniteElement>;

  template <typename FiniteElement>
  concept HascalculateLocalSystem = requires(FiniteElement fe, typename FiniteElement::FERequirementType req) {
    fe.calculateLocalSystem(req);
  };

  template <typename FiniteElement>
  concept HasFreecalculateLocalSystem = requires(FiniteElement fe, typename FiniteElement::FERequirementType req) {
    calculateLocalSystem(fe, req);
  };

  template <typename FiniteElement>
  concept MinimalFiniteElementLinearAlgebraAffordances
      = HasFreecalculateScalar<FiniteElement> || HasFreecalculateVector<FiniteElement> || HasFreecalculateMatrix<
          FiniteElement> || HascalculateScalar<FiniteElement> || HascalculateVector<FiniteElement> || HascalculateMatrix<FiniteElement>;

}  // namespace Ikarus::Concepts

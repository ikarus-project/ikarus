//
// Created by Alex on 08.04.2022.
//

#pragma once

#include <vector>
#include <set>

#include <ikarus/variables/interfaceVariable.hh>
#include <ikarus/variables/parameterFactory.hh>
namespace Ikarus {

  enum class VectorAffordances {
    noAffordance,
    forces,
    microMagneticForces
  };

  enum class MatrixAffordances {
    noAffordance,
    stiffness,
    materialstiffness,
    geometricstiffness,
    stiffnessdiffBucklingVector,
    microMagneticHessian,
    mass
  };

  enum class ScalarAffordances { noAffordance, mechanicalPotentialEnergy, microMagneticPotentialEnergy };

  struct AffordanceCollection {
    ScalarAffordances scalarAffordances{ScalarAffordances::noAffordance};
    VectorAffordances vectorAffordances{VectorAffordances::noAffordance};
    MatrixAffordances matrixAffordances{MatrixAffordances::noAffordance};
  };

  enum class ResultType { noType, magnetization, gradientNormOfMagnetization, vectorPotential, BField, HField };

  template <typename Type>
  concept FEAffordance
      = std::is_same_v<
            Type,
            ScalarAffordances> or std::is_same_v<Type, VectorAffordances> or std::is_same_v<Type, MatrixAffordances> or std::is_same_v<AffordanceCollection, MatrixAffordances>;

  std::string getResultType(const ResultType &res);

  inline constexpr VectorAffordances forces = VectorAffordances::forces;

  inline constexpr MatrixAffordances stiffness                   = MatrixAffordances::stiffness;
  inline constexpr MatrixAffordances stiffnessdiffBucklingVector = MatrixAffordances::stiffnessdiffBucklingVector;
  inline constexpr MatrixAffordances mass                        = MatrixAffordances::mass;

  inline constexpr ScalarAffordances potentialEnergy = ScalarAffordances::mechanicalPotentialEnergy;

  inline constexpr AffordanceCollection elastoStatics
      = {ScalarAffordances::mechanicalPotentialEnergy, VectorAffordances::forces, MatrixAffordances::stiffness};

  template <typename SolutionVectorType = Eigen::VectorXd, typename ParameterType = double>
  class FErequirementsBuilder;

  template <typename SolutionVectorType = Eigen::VectorXd, typename ParameterType = double>
  struct FErequirements {
    friend FErequirementsBuilder<SolutionVectorType, ParameterType>;

  public:
    const SolutionVectorType &getSolution(FESolutions &&key) const { return sols.at(key).get(); }

    const ParameterType &getParameter(FEParameter &&key) const { return parameter.at(key).get(); }

    template <FEAffordance Affordance>
    bool hasAffordance(Affordance &&affordance) const {
      if constexpr (std::is_same_v<Affordance, ScalarAffordances>)
        return affordances.scalarAffordances == affordance;
      else if constexpr (std::is_same_v<Affordance, VectorAffordances>)
        return affordances.vectorAffordances == affordance;
      else if constexpr (std::is_same_v<Affordance, MatrixAffordances>)
        return affordances.matrixAffordances == affordance;
      else if constexpr (std::is_same_v<AffordanceCollection, MatrixAffordances>)
        return affordances == affordance;
    }

  private:
    std::map<FESolutions, std::reference_wrapper<const SolutionVectorType>> sols;
    std::map<FEParameter, std::reference_wrapper<const ParameterType>> parameter;
    AffordanceCollection affordances;
  };

  template <typename SolutionVectorType, typename ParameterType>
  class FErequirementsBuilder {
  public:
    template <FEAffordance Affordance>
    FErequirementsBuilder &setAffordance(Affordance &&affordance) {
      if constexpr (std::is_same_v<Affordance, ScalarAffordances>)
        affordances.scalarAffordances = affordance;
      else if constexpr (std::is_same_v<Affordance, VectorAffordances>)
        affordances.vectorAffordances = affordance;
      else if constexpr (std::is_same_v<Affordance, MatrixAffordances>)
        affordances.matrixAffordances = affordance;
      else if constexpr (std::is_same_v<AffordanceCollection, MatrixAffordances>)
        affordances = affordance;
      return *this;
    }

    FErequirementsBuilder &setParameter(FEParameter &&key, const ParameterType &val) {
      parameter.insert({key, val});
      return *this;
    }

    FErequirementsBuilder &setSolution(FESolutions &&key, const SolutionVectorType &sol) {
      sols.insert({key, sol});
      return *this;
    }

    FErequirements<SolutionVectorType, ParameterType> build() {
      FErequirements<SolutionVectorType, ParameterType> req;
      req.sols        = sols;
      req.parameter   = parameter;
      req.affordances = affordances;
      return req;
    }

  private:
    std::map<FESolutions, std::reference_wrapper<const SolutionVectorType>> sols;
    std::map<FEParameter, std::reference_wrapper<const ParameterType>> parameter;
    AffordanceCollection affordances;
  };

  template<typename ParameterType = double>
  struct ResultTypeMap
  {
    using ResultArray = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,3,3>;
    std::map<ResultType,ResultArray> results;
  };

template <typename SolutionVectorType = Eigen::VectorXd, typename ParameterType = double>
  struct ResultRequirements : public FErequirements<SolutionVectorType,ParameterType> {

    bool isResultRequested(ResultType &&key) const { return resType.contains(key); }

    FErequirements<SolutionVectorType,ParameterType> req;
    std::set<ResultType> resType;
  };

template <typename SolutionVectorType = Eigen::VectorXd, typename ParameterType = double>
  class ResultRequirementsBuilder : public FErequirementsBuilder<SolutionVectorType,ParameterType>
  {
  public:

    template <FEAffordance Affordance>
    ResultRequirementsBuilder &setAffordance(Affordance &&affordance) {
      FErequirementsBuilder<SolutionVectorType,ParameterType>::setAffordance(affordance);
      return *this;
    }

    ResultRequirementsBuilder &setParameter(FEParameter &&key, const ParameterType &val) {
      FErequirementsBuilder<SolutionVectorType,ParameterType>::setParameter(key,val);
      return *this;
    }

    ResultRequirementsBuilder &setSolution(FESolutions &&key, const SolutionVectorType &sol) {
      FErequirementsBuilder<SolutionVectorType,ParameterType>::setSolution(key,sol);
      return *this;
    }


    ResultRequirementsBuilder & setResult(ResultType &&key) {
      resType.insert(key);
      return *this;
    }

    ResultRequirements<SolutionVectorType,ParameterType> build() {
      ResultRequirements<SolutionVectorType,ParameterType> resReq;
      resReq.req        = FErequirementsBuilder<SolutionVectorType,ParameterType>::build();
      resReq.resType   = resType;
      return resReq;
    }

  private:
    std::set<ResultType> resType;
  };

}  // namespace Ikarus
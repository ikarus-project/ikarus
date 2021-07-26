//
// Created by Alex on 27.05.2021.
//

#pragma once

#include <memory>

#include <ikarus/FiniteElements/FiniteElementPolicies.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {
  class IVariable;
  enum class VariablesTags;
}  // namespace Ikarus::Variable

namespace Ikarus::FiniteElements {
  enum class MatrixAffordances;
  enum class VectorAffordances;
  enum class ScalarAffordances;

  /** \brief A type-erased finite element */
  class IFiniteElement {
  public:
    using DofPairVectorType  = std::vector<std::pair<size_t, std::vector<Ikarus::Variable::VariablesTags>>>;
    using VariableVectorType = std::vector<Ikarus::Variable::IVariable*>;

    template <typename FE>
    explicit IFiniteElement(const FE& fe) : feimpl{std::make_unique<FEImpl<FE>>(fe)} {
      static_assert(Concepts::MinimalFiniteElementLinearAlgebraAffordances<FE>,
                    "Your element should at least provide one of the following three functions: "
                    "calculateScalar,calculateVector,calculateMatrix. These can be free or member functions.");
      static_assert(Concepts::HasSomegetEntityVariablePairs<FE>,
                    "Your element should provide the function: getEntityVariablePairs to provide degrees of freedom "
                    "definitions.");
    }

    ~IFiniteElement() = default;
    IFiniteElement(const IFiniteElement& other) : feimpl{other.feimpl->clone()} {}
    IFiniteElement& operator=(const IFiniteElement& other) {
      IFiniteElement tmp(other);
      std::swap(feimpl, tmp.feimpl);
      return *this;
    }

    IFiniteElement(IFiniteElement&&) noexcept = default;
    IFiniteElement& operator=(IFiniteElement&&) noexcept = default;

  private:
    struct FEBase {
      virtual ~FEBase()                            = default;
      virtual void do_initialize()                 = 0;
      [[nodiscard]] virtual int do_dofSize() const = 0;
      [[nodiscard]] virtual std::pair<Eigen::MatrixXd, Eigen::VectorXd> do_calculateLocalSystem(
          VariableVectorType& vars, const MatrixAffordances& matA, const VectorAffordances& vecA) const = 0;
      [[nodiscard]] virtual Eigen::MatrixXd do_calculateMatrix(VariableVectorType& vars,
                                                               const MatrixAffordances& matA) const     = 0;
      [[nodiscard]] virtual Eigen::VectorXd do_calculateVector(VariableVectorType& vars,
                                                               const VectorAffordances& vecA) const     = 0;
      [[nodiscard]] virtual double do_calculateScalar(VariableVectorType& vars,
                                                      const ScalarAffordances& scalA) const             = 0;
      [[nodiscard]] virtual DofPairVectorType do_getEntityVariablePairs() const                         = 0;
      [[nodiscard]] virtual size_t do_getEntityID() const                                               = 0;
      [[nodiscard]] virtual std::unique_ptr<FEBase> clone() const                                       = 0;
    };

    template <typename FE>
    struct FEImpl : public FEBase {
      explicit FEImpl(FE fearg) : fe{fearg} {};
      void do_initialize() final { TRYCALLFUNCTIONDONTTHROW(initialize); }
      [[nodiscard]] int do_dofSize() const final { TRYCALLFUNCTION(dofSize); }
      [[nodiscard]] std::pair<Eigen::MatrixXd, Eigen::VectorXd> do_calculateLocalSystem(
          VariableVectorType& vars, const MatrixAffordances& matA, const VectorAffordances& vecA) const final {
        TRYCALLFUNCTION(calculateLocalSystem, vars, matA, vecA);
      }
      [[nodiscard]] Eigen::MatrixXd do_calculateMatrix(VariableVectorType& vars,
                                                       const MatrixAffordances& matA) const final {
        TRYCALLFUNCTION(calculateMatrix, vars, matA);
      }
      [[nodiscard]] Eigen::VectorXd do_calculateVector(VariableVectorType& vars,
                                                       const VectorAffordances& vecA) const final {
        TRYCALLFUNCTION(calculateVector, vars, vecA);
      }
      [[nodiscard]] double do_calculateScalar(VariableVectorType& vars, const ScalarAffordances& scalA) const final {
        TRYCALLFUNCTION(calculateScalar, vars, scalA);
      }
      [[nodiscard]] DofPairVectorType do_getEntityVariablePairs() const final {
        TRYCALLFUNCTION(getEntityVariablePairs);
      }
      [[nodiscard]] size_t do_getEntityID() const final { TRYCALLFUNCTION(getEntityID); }
      [[nodiscard]] std::unique_ptr<FEBase> clone() const final { return std::make_unique<FEImpl>(*this); }
      FE fe;
    };

    std::unique_ptr<FEBase> feimpl;

    friend void initialize(IFiniteElement& fe);
    friend int dofSize(const IFiniteElement& fe);
    friend std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const IFiniteElement& fe,
                                                                            VariableVectorType& vars,
                                                                            const MatrixAffordances& matA,
                                                                            const VectorAffordances& vecA);
    friend Eigen::MatrixXd calculateMatrix(const IFiniteElement& fe, VariableVectorType& vars,
                                           const MatrixAffordances& matA);
    friend Eigen::VectorXd calculateVector(const IFiniteElement& fe, VariableVectorType& vars,
                                           const VectorAffordances& vecA);
    friend double calculateScalar(const IFiniteElement& fe, VariableVectorType& vars, const ScalarAffordances& scalA);
    friend DofPairVectorType getEntityVariablePairs(const IFiniteElement& fe);
    friend size_t getEntityID(const IFiniteElement& fe);
  };

  void initialize(IFiniteElement& fe);
  int dofSize(const IFiniteElement& fe);
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const IFiniteElement& fe,
                                                                   IFiniteElement::VariableVectorType& vars,
                                                                   const MatrixAffordances& matA,
                                                                   const VectorAffordances& vecA);
  Eigen::MatrixXd calculateMatrix(const IFiniteElement& fe, IFiniteElement::VariableVectorType& vars,
                                  const MatrixAffordances& matA);
  Eigen::MatrixXd calculateMatrix(const IFiniteElement* fe, IFiniteElement::VariableVectorType& vars,
                                  const MatrixAffordances& matA);
  Eigen::VectorXd calculateVector(const IFiniteElement& fe, IFiniteElement::VariableVectorType& vars,
                                  const VectorAffordances& vecA);
  Eigen::VectorXd calculateVector(const IFiniteElement* fe, IFiniteElement::VariableVectorType& vars,
                                  const VectorAffordances& vecA);
  double calculateScalar(const IFiniteElement& fe, IFiniteElement::VariableVectorType& vars,
                         const ScalarAffordances& scalA);
  double calculateScalar(const IFiniteElement* fe, IFiniteElement::VariableVectorType& vars,
                         const ScalarAffordances& scalA);
  IFiniteElement::DofPairVectorType getEntityVariablePairs(const IFiniteElement& fe);
  size_t getEntityID(const IFiniteElement& fe);

}  // namespace Ikarus::FiniteElements

//
// Created by Alex on 27.05.2021.
//

#pragma once
#include <memory>
#include <iostream>

#include <ikarus/FiniteElements/FiniteElementPolicies.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {
  class IVariable;
  enum class VariablesTags;
}  // namespace Ikarus::Variable
namespace Ikarus::FiniteElements {

  /** \brief A type-erased finite element */
  class IFiniteElement {
  public:
    using DofVectorType = std::vector<std::pair<size_t, std::vector<Ikarus::Variable::VariablesTags>>>;
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
      [[nodiscard]] virtual std::pair<DynMatrixd, DynVectord> do_calculateLocalSystem(
          const ElementMatrixAffordances& matA, const ElementVectorAffordances& vecA) const           = 0;
      [[nodiscard]] virtual DynMatrixd do_calculateMatrix(const ElementMatrixAffordances& matA) const = 0;
      [[nodiscard]] virtual DynVectord do_calculateVector(const ElementVectorAffordances& vecA) const = 0;
      [[nodiscard]] virtual double do_calculateScalar(const ElementScalarAffordances& scalA) const    = 0;
      [[nodiscard]] virtual DofVectorType do_getEntityVariablePairs() const                           = 0;
      [[nodiscard]] virtual std::unique_ptr<FEBase> clone() const                                     = 0;
    };

    template <typename FE>
    struct FEImpl : public FEBase {
      explicit FEImpl(FE fearg) : fe{fearg} {};
      void do_initialize() final { TRYCALLFUNCTIONDONTTHROW(initialize); }
      [[nodiscard]] int do_dofSize() const final { TRYCALLFUNCTION(dofSize); }
      [[nodiscard]] std::pair<DynMatrixd, DynVectord> do_calculateLocalSystem(
          const ElementMatrixAffordances& matA, const ElementVectorAffordances& vecA) const final {
        TRYCALLFUNCTION(calculateLocalSystem, matA, vecA);
      }
      [[nodiscard]] DynMatrixd do_calculateMatrix(const ElementMatrixAffordances& matA) const final {
        TRYCALLFUNCTION(calculateMatrix, matA);
      }
      [[nodiscard]] DynVectord do_calculateVector(const ElementVectorAffordances& vecA) const final {
        TRYCALLFUNCTION(calculateVector, vecA);
      }
      [[nodiscard]] double do_calculateScalar(const ElementScalarAffordances& scalA) const final {
        TRYCALLFUNCTION(calculateScalar, scalA);
      }
      [[nodiscard]] DofVectorType do_getEntityVariablePairs() const final { TRYCALLFUNCTION(getEntityVariablePairs); }
      [[nodiscard]] std::unique_ptr<FEBase> clone() const final { return std::make_unique<FEImpl>(*this); }
      FE fe;
    };

    std::unique_ptr<FEBase> feimpl;

    friend void initialize(IFiniteElement& fe);
    friend int dofSize(const IFiniteElement& fe);
    friend std::pair<DynMatrixd, DynVectord> calculateLocalSystem(const IFiniteElement& fe,
                                                                  const ElementMatrixAffordances& matA,
                                                                  const ElementVectorAffordances& vecA);
    friend DynMatrixd calculateMatrix(const IFiniteElement& fe, const ElementMatrixAffordances& matA);
    friend DynVectord calculateVector(const IFiniteElement& fe, const ElementVectorAffordances& vecA);
    friend double calculateScalar(const IFiniteElement& fe, const ElementScalarAffordances& scalA);
    friend DofVectorType getEntityVariablePairs(const IFiniteElement& fe);
  };

  void initialize(IFiniteElement& fe);
  int dofSize(const IFiniteElement& fe);
  std::pair<DynMatrixd, DynVectord> calculateLocalSystem(const IFiniteElement& fe, const ElementMatrixAffordances& matA,
                                                         const ElementVectorAffordances& vecA);
  DynMatrixd calculateMatrix(const IFiniteElement& fe, const ElementMatrixAffordances& matA);
  DynVectord calculateVector(const IFiniteElement& fe, const ElementVectorAffordances& vecA);
  double calculateScalar(const IFiniteElement& fe, const ElementScalarAffordances& scalA);
  IFiniteElement::DofVectorType getEntityVariablePairs(const IFiniteElement& fe);

}  // namespace Ikarus::FiniteElements

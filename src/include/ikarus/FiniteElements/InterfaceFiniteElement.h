//
// Created by Alex on 27.05.2021.
//

#pragma once
#include <ikarus/FiniteElements/FiniteElementPolicies.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {
  class IVariable;
}
namespace Ikarus::FiniteElements {

  /** \brief A type-erased finite element */
  class IFiniteElement {
  public:
    using DofVectorType = std::vector<std::pair<size_t, std::vector<Ikarus::Variable::VariablesTags>>>;
    template <typename FE>
    explicit IFiniteElement(const FE& fe) : feimpl{std::make_unique<FEImpl<FE>>(fe)} {}

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
      [[nodiscard]] virtual std::pair<DynVectord, DynMatrixd> do_calculateLocalSystem(
          const ElementVectorAffordances& vecA, const ElementMatrixAffordances& matA) const            = 0;
      [[nodiscard]] virtual DynMatrixd do_calculateMatrix(const ElementMatrixAffordances& matA) const  = 0;
      [[nodiscard]] virtual DynVectord do_calculateVector(const ElementVectorAffordances& vecA) const  = 0;
      [[nodiscard]] virtual double do_calculateScalar(const ElementScalarAffordances& scalA) const = 0;
      [[nodiscard]] virtual DofVectorType do_getEntityVariablePairs()                                  = 0;
      [[nodiscard]] virtual std::unique_ptr<FEBase> clone() const                                      = 0;
    };

    template <typename FE>
    struct FEImpl : public FEBase {
      explicit FEImpl(FE fearg) : fe{fearg} {};
      void do_initialize() final { TRYCALLFUNCTIONDONTTHROW(initialize); }
      [[nodiscard]] int do_dofSize() const final { TRYCALLFUNCTION(dofSize); }
      [[nodiscard]] std::pair<DynVectord, DynMatrixd> do_calculateLocalSystem(
          const ElementVectorAffordances& vecA, const ElementMatrixAffordances& matA) const final {
        TRYCALLFUNCTION(calculateLocalSystem, vecA, matA);
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
      [[nodiscard]] DofVectorType do_getEntityVariablePairs() final { TRYCALLFUNCTION(getEntityVariablePairs); }
      [[nodiscard]] std::unique_ptr<FEBase> clone() const final { return std::make_unique<FEImpl>(*this); }
      FE fe;
    };

    std::unique_ptr<FEBase> feimpl;

    friend void initialize(IFiniteElement& fe);
    friend int dofSize(const IFiniteElement& fe);
    friend auto calculateLocalSystem(const IFiniteElement& fe, const ElementVectorAffordances& vecA,
                                     const ElementMatrixAffordances& matA);
    friend auto calculateMatrix(const IFiniteElement& fe, const ElementMatrixAffordances& matA);
    friend auto calculateVector(const IFiniteElement& fe, const ElementVectorAffordances& vecA);
    friend auto calculateScalar(const IFiniteElement& fe, const ElementScalarAffordances& scalA);
    friend auto getEntityVariablePairs(IFiniteElement& fe);
  };

  void initialize(IFiniteElement& fe) { fe.feimpl->do_initialize(); }
  int dofSize(const IFiniteElement& fe) { return fe.feimpl->do_dofSize(); }
  auto calculateLocalSystem(const IFiniteElement& fe, const ElementVectorAffordances& vecA,
                            const ElementMatrixAffordances& matA) {
    return fe.feimpl->do_calculateLocalSystem(vecA, matA);
  }
  auto calculateMatrix(const IFiniteElement& fe, const ElementMatrixAffordances& matA) {
    return fe.feimpl->do_calculateMatrix(matA);
  }
  auto calculateVector(const IFiniteElement& fe, const ElementVectorAffordances& vecA) {
    return fe.feimpl->do_calculateVector(vecA);
  }
  auto calculateScalar(const IFiniteElement& fe, const ElementScalarAffordances& scalA) {
    return fe.feimpl->do_calculateScalar(scalA);
  }
  auto getEntityVariablePairs(IFiniteElement& fe) { return fe.feimpl->do_getEntityVariablePairs(); }

}  // namespace Ikarus::FiniteElements

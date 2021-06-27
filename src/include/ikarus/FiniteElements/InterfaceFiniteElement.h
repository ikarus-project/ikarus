//
// Created by Alex on 27.05.2021.
//

#pragma once
#include <ikarus/FiniteElements/PhysicalElementPolicies.h>

namespace Ikarus::Variable {
  class GenericVariable;
}
namespace Ikarus::FiniteElements {

  /** \brief A type-erased finite element */
  class IFiniteElement {
  public:
    using DofVectorType = std::vector<Ikarus::Variable::GenericVariable*>;
    template <typename FE>
    explicit IFiniteElement(const FE& fe) : feimpl{std::make_unique<FEImpl<FE> >(fe)} {}

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
      virtual ~FEBase()                                                                       = default;
      virtual void do_initialize()                                                            = 0;
      [[nodiscard]] virtual int do_dofSize() const                                            = 0;
      [[nodiscard]] virtual std::pair<DynVectord, DynMatrixd> do_calculateLocalSystem() const = 0;
      [[nodiscard]] virtual DynMatrixd do_calculateLHS() const                                = 0;
      [[nodiscard]] virtual DynVectord do_calculateRHS() const                                = 0;
      [[nodiscard]] virtual DofVectorType do_getDofVector()                                   = 0;
      [[nodiscard]] virtual std::unique_ptr<FEBase> clone() const = 0;
    };

    template <typename FE>
    struct FEImpl : public FEBase {
      explicit FEImpl(FE fearg) : fe{fearg} {};
      void do_initialize() final { TRYCALLFUNCTIONDONTTHROW(initialize); }
      [[nodiscard]] int do_dofSize() const final { TRYCALLFUNCTION(dofSize); }
      [[nodiscard]] std::pair<DynVectord, DynMatrixd> do_calculateLocalSystem() const final {
        TRYCALLFUNCTION(calculateLocalSystem);
      }
      [[nodiscard]] DynMatrixd do_calculateLHS() const final { TRYCALLFUNCTION(calculateLHS); }
      [[nodiscard]] DynVectord do_calculateRHS() const final { TRYCALLFUNCTION(calculateRHS); }
      [[nodiscard]] DofVectorType do_getDofVector() final { TRYCALLFUNCTION(getDofVector); }
      [[nodiscard]] std::unique_ptr<FEBase> clone() const final { return std::make_unique<FEImpl>(*this); }
      FE fe;
    };

    std::unique_ptr<FEBase> feimpl;

    friend void initialize(IFiniteElement& fe);
    friend int dofSize(const IFiniteElement& fe);
    friend auto calculateLocalSystem(const IFiniteElement& fe);
    friend auto calculateLHS(const IFiniteElement& fe);
    friend auto calculateRHS(const IFiniteElement& fe);
    friend auto getDofVector(IFiniteElement& fe);
  };

  void initialize(IFiniteElement& fe) { fe.feimpl->do_initialize(); }
  int dofSize(const IFiniteElement& fe) { return fe.feimpl->do_dofSize(); }
  auto calculateLocalSystem(const IFiniteElement& fe) { return fe.feimpl->do_calculateLocalSystem(); }
  auto calculateLHS(const IFiniteElement& fe) { return fe.feimpl->do_calculateLHS(); }
  auto calculateRHS(const IFiniteElement& fe) { return fe.feimpl->do_calculateRHS(); }
  auto getDofVector(IFiniteElement& fe) { return fe.feimpl->do_getDofVector(); }

}  // namespace Ikarus::FiniteElements

//
// Created by Alex on 27.05.2021.
//

#pragma once

#include <ikarus/FiniteElements/FiniteElementInterface.h>
#include <ikarus/FiniteElements/PhysicalElementPolicies.h>

namespace Ikarus::PhysicalElements {

  /** \brief A type-erased physical element */
  class GenericFE {
  public:
    template <typename FE> explicit GenericFE(const FE& fe) : feimpl{std::make_unique<FEImpl<FE> >(fe)} {}

    ~GenericFE() = default;
    GenericFE(const GenericFE& other) : feimpl{other.feimpl->clone()} {}
    GenericFE& operator=(const GenericFE& other) {
      GenericFE tmp(other);  // Temporary-swap idiom
      std::swap(feimpl, tmp.feimpl);
      return *this;
    }

    GenericFE(GenericFE&&) noexcept = default;
    GenericFE& operator=(GenericFE&&) noexcept = default;

  private:
    struct FEBase {
      virtual ~FEBase() = default;
      virtual void do_initialize() = 0;
      [[nodiscard]] virtual int do_dofSize() const = 0;
      [[nodiscard]] virtual std::pair<DynVectord, DynMatrixd> do_calculateLocalSystem() const = 0;
      [[nodiscard]] virtual DynMatrixd do_calculateLHS() const = 0;
      [[nodiscard]] virtual DynVectord do_calculateRHS() const = 0;
      [[nodiscard]] virtual DynArrayXi do_getDofVector() const = 0;
      [[nodiscard]] virtual std::unique_ptr<FEBase> clone() const = 0;  // Prototype Design Pattern
    };

    template <typename FE> struct FEImpl : public FEBase {
      explicit FEImpl(FE fearg) : fe{fearg} {};
      void do_initialize() final { TRYCALLMEMBERFUNCTIONDONTTHROW(initialize); }
      [[nodiscard]] int do_dofSize() const final { TRYCALLMEMBERFUNCTION(dofSize); }
      [[nodiscard]] std::pair<DynVectord, DynMatrixd> do_calculateLocalSystem() const final {
        TRYCALLMEMBERFUNCTION(calculateLocalSystem);
      }
      [[nodiscard]] DynMatrixd do_calculateLHS() const final { TRYCALLMEMBERFUNCTION(calculateLHS); }
      [[nodiscard]] DynVectord do_calculateRHS() const final { TRYCALLMEMBERFUNCTION(calculateRHS); };
      [[nodiscard]] DynArrayXi do_getDofVector() const final { TRYCALLMEMBERFUNCTION(getDofVector); };
      [[nodiscard]] std::unique_ptr<FEBase> clone() const final { return std::make_unique<FEImpl>(*this); }
      FE fe;
    };

    std::unique_ptr<FEBase> feimpl;  // Pimpl idiom / Bridge Design Patterns

    friend void initialize(GenericFE& fe);
    friend int dofSize(const GenericFE& fe);
    friend auto calculateLocalSystem(const GenericFE& fe);
    friend auto calculateLHS(const GenericFE& fe);
    friend auto calculateRHS(const GenericFE& fe);
    friend auto getDofVector(const GenericFE& fe);
  };

  void initialize(GenericFE& fe) { fe.feimpl->do_initialize(); }
  int dofSize(const GenericFE& fe) { return fe.feimpl->do_dofSize(); }
  auto calculateLocalSystem(const GenericFE& fe) { return fe.feimpl->do_calculateLocalSystem(); }
  auto calculateLHS(const GenericFE& fe) { return fe.feimpl->do_calculateLHS(); }
  auto calculateRHS(const GenericFE& fe) { return fe.feimpl->do_calculateRHS(); }
  auto getDofVector(const GenericFE& fe) { return fe.feimpl->do_getDofVector(); }

}  // namespace Ikarus::PhysicalElements

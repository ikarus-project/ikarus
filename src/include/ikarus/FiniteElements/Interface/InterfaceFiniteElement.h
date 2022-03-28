//
// Created by Alex on 27.05.2021.
//

#pragma once

#include <memory>
#include <optional>

#include "FiniteElementFunctionConcepts.h"
#include "ikarus/utils/LinearAlgebraTypedefs.h"

namespace Ikarus::Variable {
  class IVariable;
  enum class VariableTags;
}  // namespace Ikarus::Variable

namespace Ikarus::FiniteElements {
  enum class MatrixAffordances;
  enum class VectorAffordances;
  enum class ScalarAffordances;

  /** \brief A type-erased finite element but templated over the localView and the passed SolutionType*/
  template <typename LocalView, typename SolutionVectorType>
  class IFiniteElement {
  public:
    using DofPairVectorType = std::vector<VariableIndicesPair>;
    using FERequirementType = FErequirements<SolutionVectorType>;
    using GlobalIndex       = typename LocalView::MultiIndex;
    //    using DataVectorType     = typename std::optional<std::reference_wrapper<VariableVectorType>>;

    template <typename FE>
    explicit IFiniteElement(const FE &fe) : feimpl{std::make_unique<FEImpl<FE>>(fe)} {
      static_assert(Concepts::MinimalFiniteElementLinearAlgebraAffordances<FE>,
                    "Your element should at least provide one of the following three functions: "
                    "calculateScalar,calculateVector,calculateMatrix. These can be free or member functions.");
    }

    ~IFiniteElement() = default;
    IFiniteElement(const IFiniteElement &other) : feimpl{other.feimpl->clone()} {}
    IFiniteElement &operator=(const IFiniteElement &other) {
      IFiniteElement tmp(other);
      std::swap(feimpl, tmp.feimpl);
      return *this;
    }

    IFiniteElement(IFiniteElement &&) noexcept = default;
    IFiniteElement &operator=(IFiniteElement &&) noexcept = default;

  private:
    struct FEBase {
      virtual ~FEBase() = default;
      [[nodiscard]] virtual std::pair<Eigen::MatrixXd, Eigen::VectorXd> do_calculateLocalSystem(
          const FERequirementType &par) const                                                      = 0;
      [[nodiscard]] virtual Eigen::MatrixXd do_calculateMatrix(const FERequirementType &par) const = 0;
      //      virtual void do_bind(const LocalViewEmbedded &localView) const = 0;
      [[nodiscard]] virtual Eigen::VectorXd do_calculateVector(const FERequirementType &par) const = 0;
      [[nodiscard]] virtual double do_calculateScalar(const FERequirementType &par) const          = 0;
      [[nodiscard]] virtual std::vector<GlobalIndex> do_globalIndices() const                      = 0;
      [[nodiscard]] virtual std::unique_ptr<FEBase> clone() const                                  = 0;
    };

    template <typename FE>
    struct FEImpl : public FEBase {
      explicit FEImpl(FE fearg) : fe{fearg} {};
      [[nodiscard]] std::pair<Eigen::MatrixXd, Eigen::VectorXd> do_calculateLocalSystem(
          const FERequirementType &par) const final {
        TRYCALLFUNCTIONANDRETURN(calculateLocalSystem, par);
      }
      [[nodiscard]] Eigen::MatrixXd do_calculateMatrix(const FERequirementType &par) const final {
        TRYCALLFUNCTIONANDRETURN(calculateMatrix, par);
      }
      [[nodiscard]] Eigen::VectorXd do_calculateVector(const FERequirementType &par) const final {
        TRYCALLFUNCTIONANDRETURN(calculateVector, par);
      }
      [[nodiscard]] double do_calculateScalar(const FERequirementType &par) const final {
        TRYCALLFUNCTIONANDRETURN(calculateScalar, par);
      }
      [[nodiscard]] std::vector<GlobalIndex> do_globalIndices() const final { TRYCALLFUNCTIONANDRETURN(globalIndices); }

      [[nodiscard]] std::unique_ptr<FEBase> clone() const final { return std::make_unique<FEImpl>(*this); }
      FE fe;
    };

    std::unique_ptr<FEBase> feimpl;

  public:
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const FERequirementType &par) const;
    Eigen::MatrixXd calculateMatrix(const FERequirementType &par) const;
    Eigen::VectorXd calculateVector(const FERequirementType &par) const;
    double calculateScalar(const FERequirementType &par) const;
    std::vector<GlobalIndex> globalIndices() const;
  };

#include "InterfaceFiniteElement.inl"

}  // namespace Ikarus::FiniteElements

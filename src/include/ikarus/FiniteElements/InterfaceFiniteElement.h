//
// Created by Alex on 27.05.2021.
//

#pragma once

#include <memory>
#include <optional>

#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/Grids/EntityHelperFunctions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {
  class IVariable;
  enum class VariableTags;
}  // namespace Ikarus::Variable

namespace Ikarus::FiniteElements {
  enum class MatrixAffordances;
  enum class VectorAffordances;
  enum class ScalarAffordances;

  /** \brief A type-erased finite element */
  class IFiniteElement {
  public:
    using DofPairVectorType = std::vector<VariableIndicesPair>;
    using FERequirementType = FErequirements;
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
      virtual ~FEBase()                            = default;
      [[nodiscard]] virtual int do_dofSize() const = 0;
      [[nodiscard]] virtual std::pair<Eigen::MatrixXd, Eigen::VectorXd> do_calculateLocalSystem(
          const FERequirementType &par) const                                                      = 0;
      [[nodiscard]] virtual Eigen::MatrixXd do_calculateMatrix(const FERequirementType &par) const = 0;
      [[nodiscard]] virtual Eigen::VectorXd do_calculateVector(const FERequirementType &par) const = 0;
      [[nodiscard]] virtual double do_calculateScalar(const FERequirementType &par) const          = 0;
      [[nodiscard]] virtual DofPairVectorType do_getEntityVariableTuple() const                    = 0;
      [[nodiscard]] virtual unsigned int do_subEntities(unsigned int codim) const                  = 0;
      [[nodiscard]] virtual size_t do_subIndex(int i, unsigned int codim) const                    = 0;
      [[nodiscard]] virtual unsigned int do_dimension() const                                      = 0;
      [[nodiscard]] virtual std::unique_ptr<FEBase> clone() const                                  = 0;
    };

    template <typename FE>
    struct FEImpl : public FEBase {
      explicit FEImpl(FE fearg) : fe{fearg} {};
      [[nodiscard]] int do_dofSize() const final { TRYCALLFUNCTION(dofSize); }
      [[nodiscard]] std::pair<Eigen::MatrixXd, Eigen::VectorXd> do_calculateLocalSystem(
          const FERequirementType &par) const final {
        TRYCALLFUNCTION(calculateLocalSystem, par);
      }
      [[nodiscard]] Eigen::MatrixXd do_calculateMatrix(const FERequirementType &par) const final {
        TRYCALLFUNCTION(calculateMatrix, par);
      }
      [[nodiscard]] Eigen::VectorXd do_calculateVector(const FERequirementType &par) const final {
        TRYCALLFUNCTION(calculateVector, par);
      }
      [[nodiscard]] double do_calculateScalar(const FERequirementType &par) const final {
        TRYCALLFUNCTION(calculateScalar, par);
      }
      [[nodiscard]] DofPairVectorType do_getEntityVariableTuple() const final {
        TRYCALLFUNCTION(getEntityVariableTuple);
      }
      [[nodiscard]] unsigned int do_subEntities(unsigned int codim) const final { TRYCALLFUNCTION(subEntities, codim); }
      [[nodiscard]] size_t do_subIndex(int i, unsigned int codim) const final { TRYCALLFUNCTION(subIndex, i, codim); }
      [[nodiscard]] unsigned int do_dimension() const final { TRYCALLFUNCTION(dimension); }
      [[nodiscard]] std::unique_ptr<FEBase> clone() const final { return std::make_unique<FEImpl>(*this); }
      FE fe;
    };

    std::unique_ptr<FEBase> feimpl;

    friend void initialize(IFiniteElement &fe);
    friend int dofSize(const IFiniteElement &fe);
    friend std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const IFiniteElement &fe,
                                                                            const FERequirementType &par);

    friend Eigen::MatrixXd calculateMatrix(const IFiniteElement &fe, const FERequirementType &par);
    friend Eigen::VectorXd calculateVector(const IFiniteElement &fe, const FERequirementType &par);
    friend double calculateScalar(const IFiniteElement &fe, const FERequirementType &par);
    friend DofPairVectorType getEntityVariableTuple(const IFiniteElement &fe);
    friend unsigned int subEntities(const IFiniteElement &fe, unsigned int codim);
    friend unsigned int dimension(const IFiniteElement &fe);
    friend size_t subIndex(const IFiniteElement &fe, int i, unsigned int codim);
  };

  int dofSize(const IFiniteElement &fe);
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const IFiniteElement &fe,
                                                                   const IFiniteElement::FERequirementType &par);
  Eigen::MatrixXd calculateMatrix(const IFiniteElement &fe, const IFiniteElement::FERequirementType &part);
  Eigen::VectorXd calculateVector(const IFiniteElement &fe, const IFiniteElement::FERequirementType &par);
  double calculateScalar(const IFiniteElement &fe, const IFiniteElement::FERequirementType &par);
  IFiniteElement::DofPairVectorType getEntityVariableTuple(const IFiniteElement &fe);
  unsigned int subEntities(const IFiniteElement &fe, unsigned int codim);
  size_t subIndex(const IFiniteElement &fe, int i, unsigned int codim);
  unsigned int dimension(const IFiniteElement &fe);

}  // namespace Ikarus::FiniteElements

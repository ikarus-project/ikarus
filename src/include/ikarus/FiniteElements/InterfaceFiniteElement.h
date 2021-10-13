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
  enum class VariablesTags;
}  // namespace Ikarus::Variable

namespace Ikarus::FiniteElements {
  enum class MatrixAffordances;
  enum class VectorAffordances;
  enum class ScalarAffordances;

  /** \brief A type-erased finite element */
  class IFiniteElement {
  public:
    using DofPairVectorType  = std::vector<DofAtEntity>;
    using VariableVectorType = FEValues;
    using DataVectorType     = typename std::optional<std::reference_wrapper<VariableVectorType>>;

    template <typename FE>
    explicit IFiniteElement(const FE &fe) : feimpl{std::make_unique<FEImpl<FE>>(fe)} {
      static_assert(Concepts::MinimalFiniteElementLinearAlgebraAffordances<FE>,
                    "Your element should at least provide one of the following three functions: "
                    "calculateScalar,calculateVector,calculateMatrix. These can be free or member functions.");
      //    static_assert(Concepts::HasSomegetEntityVariableTuple<FE>,
      //                  "Your element should provide the function: getEntityVariableTuple to provide degrees of
      //                  freedom " "definitions.");
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
      virtual void do_initialize()                 = 0;
      [[nodiscard]] virtual int do_dofSize() const = 0;
      [[nodiscard]] virtual std::pair<Eigen::MatrixXd, Eigen::VectorXd> do_calculateLocalSystem(
          const MatrixAffordances &matA, const VectorAffordances &vecA, VariableVectorType &vars,
          DataVectorType &) const                                                          = 0;
      [[nodiscard]] virtual Eigen::MatrixXd do_calculateMatrix(const MatrixAffordances &matA, VariableVectorType &vars,
                                                               DataVectorType &data) const = 0;
      [[nodiscard]] virtual Eigen::VectorXd do_calculateVector(const VectorAffordances &vecA, VariableVectorType &vars,
                                                               DataVectorType &data) const = 0;
      [[nodiscard]] virtual double do_calculateScalar(const ScalarAffordances &scalA, VariableVectorType &vars,
                                                      DataVectorType &data) const          = 0;
      [[nodiscard]] virtual DofPairVectorType do_getEntityVariableTuple() const            = 0;
      [[nodiscard]] virtual unsigned int do_subEntities(unsigned int codim) const          = 0;
      [[nodiscard]] virtual size_t do_subIndex(int i, unsigned int codim) const            = 0;
      [[nodiscard]] virtual unsigned int do_dimension() const                              = 0;
      [[nodiscard]] virtual std::unique_ptr<FEBase> clone() const                          = 0;
    };

    template <typename FE>
    struct FEImpl : public FEBase {
      explicit FEImpl(FE fearg) : fe{fearg} {};
      void do_initialize() final { TRYCALLFUNCTIONDONTTHROW(initialize); }
      [[nodiscard]] int do_dofSize() const final { TRYCALLFUNCTION(dofSize); }
      [[nodiscard]] std::pair<Eigen::MatrixXd, Eigen::VectorXd> do_calculateLocalSystem(const MatrixAffordances &matA,
                                                                                        const VectorAffordances &vecA,
                                                                                        VariableVectorType &vars,
                                                                                        DataVectorType &data
                                                                                        = std::nullopt) const final {
        if (data) {
          TRYCALLFUNCTION(calculateLocalSystem, matA, vecA, vars, data);
        } else
          TRYCALLFUNCTIONWITHOUTDATA(calculateLocalSystem, matA, vecA, vars);
      }
      [[nodiscard]] Eigen::MatrixXd do_calculateMatrix(const MatrixAffordances &matA, VariableVectorType &vars,
                                                       DataVectorType &data = std::nullopt) const final {
        if (data) {
          TRYCALLFUNCTION(calculateMatrix, matA, vars, data);
        } else
          TRYCALLFUNCTIONWITHOUTDATA(calculateMatrix, matA, vars);
      }
      [[nodiscard]] Eigen::VectorXd do_calculateVector(const VectorAffordances &vecA, VariableVectorType &vars,
                                                       DataVectorType &data = std::nullopt) const final {
        if (data) {
          TRYCALLFUNCTION(calculateVector, vecA, vars, data);
        } else
          TRYCALLFUNCTIONWITHOUTDATA(calculateVector, vecA, vars);
      }
      [[nodiscard]] double do_calculateScalar(const ScalarAffordances &scalA, VariableVectorType &vars,
                                              DataVectorType &data = std::nullopt) const final {
        if (data) {
          TRYCALLFUNCTION(calculateScalar, scalA, vars, data);
        } else
          TRYCALLFUNCTIONWITHOUTDATA(calculateScalar, scalA, vars);
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
                                                                            const MatrixAffordances &matA,
                                                                            const VectorAffordances &vecA,
                                                                            VariableVectorType &vars,
                                                                            DataVectorType data);

    friend Eigen::MatrixXd calculateMatrix(const IFiniteElement &fe, const MatrixAffordances &matA,
                                           VariableVectorType &vars, DataVectorType data);
    friend Eigen::VectorXd calculateVector(const IFiniteElement &fe, const VectorAffordances &vecA,
                                           VariableVectorType &vars, DataVectorType data);
    friend double calculateScalar(const IFiniteElement &fe, const ScalarAffordances &scalA, VariableVectorType &vars,
                                  DataVectorType data);
    friend DofPairVectorType getEntityVariableTuple(const IFiniteElement &fe);
    friend unsigned int subEntities(const IFiniteElement &fe, unsigned int codim);
    friend unsigned int dimension(const IFiniteElement &fe);
    friend size_t subIndex(const IFiniteElement &fe, int i, unsigned int codim);
  };

  void initialize(IFiniteElement &fe);
  int dofSize(const IFiniteElement &fe);
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const IFiniteElement &fe,
                                                                   const MatrixAffordances &matA,
                                                                   const VectorAffordances &vecA,
                                                                   IFiniteElement::VariableVectorType &vars,
                                                                   IFiniteElement::DataVectorType data = std::nullopt);
  Eigen::MatrixXd calculateMatrix(const IFiniteElement &fe, const MatrixAffordances &matA,
                                  IFiniteElement::VariableVectorType &vars,
                                  IFiniteElement::DataVectorType data = std::nullopt);
  Eigen::VectorXd calculateVector(const IFiniteElement &fe, const VectorAffordances &vecA,
                                  IFiniteElement::VariableVectorType &vars,
                                  IFiniteElement::DataVectorType data = std::nullopt);
  double calculateScalar(const IFiniteElement &fe, const ScalarAffordances &scalA,
                         IFiniteElement::VariableVectorType &vars, IFiniteElement::DataVectorType data = std::nullopt);
  IFiniteElement::DofPairVectorType getEntityVariableTuple(const IFiniteElement &fe);
  unsigned int subEntities(const IFiniteElement &fe, unsigned int codim);
  size_t subIndex(const IFiniteElement &fe, int i, unsigned int codim);
  unsigned int dimension(const IFiniteElement &fe);

}  // namespace Ikarus::FiniteElements

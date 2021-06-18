//
// Created by Alex on 25.05.2021.
//

#pragma once

#include <span>
#include <memory>
#include <Eigen/Core>

#include "VariableInterface.h"

#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Variable {

    class GenericVariable {
    public:
        template<Concepts::Variable VAR>
        explicit GenericVariable(const VAR &vo)
                :variableImpl{std::make_unique < VarImpl < VAR > > (vo)} {}

        ~GenericVariable() = default;

        GenericVariable(const GenericVariable &other)
                : variableImpl{other.variableImpl->clone()} {}

        GenericVariable &operator=(const GenericVariable &other) {
            GenericVariable tmp(other);  // Temporary-swap idiom
            std::swap(variableImpl, tmp.variableImpl);
            return *this;
        }

        GenericVariable(GenericVariable &&) noexcept = default;

        GenericVariable &operator=(GenericVariable &&) noexcept = default;

    private:
        struct VarBase {
            virtual ~VarBase() = default;

            [[nodiscard]] virtual int do_valueSize() const = 0;
            [[nodiscard]] virtual int do_correctionSize() const = 0;
            [[nodiscard]] virtual bool do_equalComparison(const GenericVariable &other) const = 0;
            [[nodiscard]] virtual bool do_lessComparison(const GenericVariable &other) const = 0;
            virtual void do_update(const Eigen::Ref<const Ikarus::DynVectord> &other) = 0;
            virtual void do_setValue(const Eigen::Ref<const Ikarus::DynVectord> &other) = 0;
            [[nodiscard]] virtual Ikarus::DynVectord do_getValue() const = 0;
            [[nodiscard]] virtual size_t do_getTag() const = 0;
            [[nodiscard]] virtual std::unique_ptr <VarBase> clone() const = 0; // Prototype Design Pattern
        };

        template<typename VAR>
        struct VarImpl : public VarBase {
            explicit VarImpl(VAR voarg) : vo{voarg} {};

            [[nodiscard]]  int do_valueSize() const final { return VAR::valueSize; }
            [[nodiscard]]  int do_correctionSize() const final { return VAR::correctionSize; }
            [[nodiscard]] size_t do_getTag() const final { return vo.getTag(); };
            [[nodiscard]] Ikarus::DynVectord do_getValue() const final { return vo.getValue(); };
            void do_update(const Eigen::Ref<const Ikarus::DynVectord> &other) final { return vo.update(other); }
            void do_setValue(const Eigen::Ref<const Ikarus::DynVectord> &other) final { return vo.setValue(other); }
            [[nodiscard]] bool do_equalComparison(const GenericVariable &other) const final {
                return (this->do_getTag() == getTag(other));
            };

            [[nodiscard]] bool do_lessComparison(const GenericVariable &other) const final {
                return (this->do_getTag() < getTag(other));
            };

            [[nodiscard]] std::unique_ptr <VarBase> clone() const final { return std::make_unique<VarImpl>(*this); }

            VAR vo;
        };

        std::unique_ptr <VarBase> variableImpl; // Pimpl idiom / Bridge Design Pattern

        friend void update(GenericVariable &vo, const Eigen::Ref<const Ikarus::DynVectord> &correction);
        friend void setValue(GenericVariable &vo, const Eigen::Ref<const Ikarus::DynVectord> &value);
        friend Ikarus::DynVectord getValue(const GenericVariable &vo);
        friend int valueSize(const GenericVariable &vo);
        friend int correctionSize(const GenericVariable &vo);
        friend bool operator==(const GenericVariable &var, const GenericVariable &other);
        friend bool operator<(const GenericVariable &var, const GenericVariable &other);
        friend size_t getTag(const GenericVariable &var);
        friend std::ostream &operator<<(std::ostream &s, const GenericVariable &var);
    };

    void update(GenericVariable &vo, const Eigen::Ref<const Ikarus::DynVectord> &correction);
    void setValue(GenericVariable &vo, const Eigen::Ref<const Ikarus::DynVectord> &value);
    Ikarus::DynVectord getValue(const GenericVariable &vo);
    int valueSize(const GenericVariable &vo);
    int correctionSize(const GenericVariable &vo);
    bool operator==(const GenericVariable &var, const GenericVariable &other);
    bool operator<(const GenericVariable &var, const GenericVariable &other);
    size_t getTag(const GenericVariable &var);
    std::ostream &operator<<(std::ostream &s, const GenericVariable &var);
    size_t valueSize(std::span<const GenericVariable> varSpan);
    size_t correctionSize(std::span<const GenericVariable> varSpan);
    void update(std::span<GenericVariable> varSpan,const Ikarus::DynVectord& correction);
}

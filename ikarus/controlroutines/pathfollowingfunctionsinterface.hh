// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file pathfollowingfunctionsinterface.hh
 * \brief

 */
 #pragma once 
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
 namespace Ikarus{
template<typename NLO>
class SubsidaryFunction {

    struct SFConcept {
        virtual ~SFConcept() = default;
        virtual std::unique_ptr<SFConcept> copy_() const = 0;
        virtual void call_(SubsidiaryArgs& args) const = 0;
        virtual void initialPrediction_(NLO& nonLinearOperator, SubsidiaryArgs& args) = 0;
        virtual void intermediatePrediction_(NLO& nonLinearOperator, SubsidiaryArgs& args) = 0;
        virtual std::string name_() const = 0;
    };

    template <typename T> struct SFModel : SFConcept {
        SFModel(const T& t) : value_(t) {}
        std::unique_ptr<SFConcept> copy_() const override {
            return std::make_unique<SFModel>(value_);
        }
        void call_(SubsidiaryArgs& args) const override { value_.call(args); }
        void initialPrediction_(NLO& nonLinearOperator, SubsidiaryArgs& args) override { value_.initialPrediction(nonLinearOperator, args); }
        void intermediatePrediction_(NLO& nonLinearOperator, SubsidiaryArgs& args) override { value_.intermediatePrediction(nonLinearOperator, args); }
        std::string name_() const override { return value_.name(); }

        T value_;
    };

public:
    template <typename T> SubsidaryFunction(const T &t) {
        value_ = std::make_unique<SFModel<T>>(t);
    }
    SubsidaryFunction(const SubsidaryFunction& other) : value_(other.value_->copy_()) {}

    void operator()(SubsidiaryArgs& args) const { value_->call_(args); }
    void initialPrediction(NLO& nonLinearOperator, SubsidiaryArgs& args) { value_->initialPrediction_(nonLinearOperator, args); }
    void intermediatePrediction(NLO& nonLinearOperator, SubsidiaryArgs& args) { value_->intermediatePrediction_(nonLinearOperator, args); }
    std::string name() const { return value_->name_(); }

private:
    std::unique_ptr<SFConcept> value_;
};

}
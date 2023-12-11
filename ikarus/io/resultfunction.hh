// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <type_traits>
#include <utility>

#include <dune/vtk/vtkwriter.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/io/resultevaluators.hh>
#include <ikarus/utils/eigendunetransformations.hh>

namespace Ikarus {
  namespace Impl {
    struct DefaultUserFunction {};
  }  // namespace Impl
  /**
   * \brief Wrapper to evaluate results for a vtkwriter.
   * Usage:
   *   auto resReq = Ikarus::ResultRequirements()
                                .insertGlobalSolution(Ikarus::FESolutions::displacement, d)
                                .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                                .addResultRequest(ResultType::PK2Stress);
    auto resultFunction = std::make_shared<ResultFunction<ElementType>>(&fes, resReq);
   * Usage with dune-vtk : vtkWriter2.addPointData(Dune::Vtk::Function<GridView>( resultFunction));
   * or with Dune vtkWriter.addVertexData(resultFunction);
   */

  template <typename ElementType_, typename UserFunction = Impl::DefaultUserFunction>
  class ResultFunction : public Dune::VTKFunction<typename ElementType_::GridView> {
  public:
    using ElementType            = ElementType_;
    using ResultRequirements     = typename ElementType::ResultRequirementsType;
    using GridView               = typename ElementType::GridView;
    using ctype                  = typename GridView::ctype;
    constexpr static int griddim = GridView::dimension;
    typedef typename GridView::template Codim<0>::Entity Entity;

    double evaluate(int comp, const Entity& e, const Dune::FieldVector<ctype, griddim>& local) const override {
      auto index = gridView.indexSet().index(e);
      return evaluateStressComponent(index, local, comp);
    }

    [[nodiscard]] int ncomps() const override {
      if constexpr (std::is_same_v<UserFunction, Impl::DefaultUserFunction>) {
        Dune::FieldVector<ctype, griddim> val(0.0);

        fes_->at(0).calculateAt(resultRequirements_, val, resultTypeMap);
        if (resultRequirements_.getRequestedResult() != resultTypeMap.getSingleResult().first)
          DUNE_THROW(Dune::InvalidStateException, "The return result should be the requested one");

        auto sigma = resultTypeMap.getSingleResult().second;

        return static_cast<int>(sigma.rows() * sigma.cols());
      } else
        return userFunction_.ncomps();
    }

    [[nodiscard]] constexpr std::string name() const override {
      if constexpr (std::is_same_v<UserFunction, Impl::DefaultUserFunction>)
        return toString(resultRequirements_.getRequestedResult());
      else
        return userFunction_.name();
    }

    ResultFunction(std::vector<ElementType>* fes, const ResultRequirements& req, UserFunction userFunction = {})
        : gridView{fes->at(0).localView().globalBasis().gridView()},
          resultRequirements_{req},
          fes_{fes},
          userFunction_{userFunction} {
      if constexpr (!std::is_same_v<UserFunction, Impl::DefaultUserFunction>) userFunction_ = userFunction;
    }

  private:
    double evaluateStressComponent(int eleID, const Dune::FieldVector<ctype, griddim>& local, int comp) const {
      if constexpr (!std::is_same_v<UserFunction, Impl::DefaultUserFunction>)
        return userFunction_(fes_->at(eleID), resultRequirements_, local, comp);
      else {
        fes_->at(eleID).calculateAt(resultRequirements_, local, resultTypeMap);
        auto result = resultTypeMap.getSingleResult().second;
        return result(comp);
      }
    }

    GridView gridView;
    ResultRequirements resultRequirements_;
    std::vector<ElementType>* fes_;
    mutable ResultTypeMap<ctype> resultTypeMap;
    [[no_unique_address]] std::string name_{};
    UserFunction userFunction_;
  };
}  // namespace Ikarus

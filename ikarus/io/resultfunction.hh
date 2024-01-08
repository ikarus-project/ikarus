// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file resultfunction.hh
 * @brief Ikarus Result Evaluators for Stress Analysis
 * @ingroup io
 *
 */

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
   * @brief Wrapper to evaluate results for a vtkwriter.
   * @details
   * Usage:
   * @code
   *   auto resReq = Ikarus::ResultRequirements()
                                .insertGlobalSolution(Ikarus::FESolutions::displacement, d)
                                .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                                .addResultRequest(ResultType::PK2Stress);
    auto resultFunction = std::make_shared<ResultFunction<ElementType>>(&fes, resReq);

   * vtkWriter.addPointData(Dune::Vtk::Function<GridView>( resultFunction));
   * // or with Dunes native Vtk
   * vtkWriter.addVertexData(resultFunction);
     * @endcode
  * @ingroup io
  * @tparam ElementType_ Type of the finite element
   * @tparam UserFunction Type of the user-defined function for custom result evaluation (default is
  DefaultUserFunction)
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

    /**
     * @brief Evaluate the component at a given entity and local coordinates.
     *
     * This function is required by the Dune::VTKFunction interface.
     *
     * @param comp Stress component index
     * @param e Entity on which to evaluate the stress
     * @param local Local coordinates within the entity
     * @return Stress component value
     */
    double evaluate(int comp, const Entity& e, const Dune::FieldVector<ctype, griddim>& local) const override {
      auto index = gridView.indexSet().index(e);
      return evaluateComponent(index, local, comp);
    }

    /**
     * @brief Get the number of components.
     *
     * This function is required by the Dune::VTKFunction interface.
     *
     * @return Number of stress components
     */
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

    /**
     * @brief Get the name of the result type.
     *
     * This function is required by the Dune::VTKFunction interface.
     *
     * @return String representing the name of the result type
     */
    [[nodiscard]] constexpr std::string name() const override {
      if constexpr (std::is_same_v<UserFunction, Impl::DefaultUserFunction>)
        return toString(resultRequirements_.getRequestedResult());
      else
        return userFunction_.name();
    }

    /**
     * @brief Constructor for ResultFunction.
     *
     * Constructs a ResultFunction object with given finite elements, result requirements, and an optional user
     * function.
     *
     * @param fes Pointer to a vector of finite elements
     * @param req Result requirements for evaluation
     * @param userFunction User-defined function for custom result evaluation (default is DefaultUserFunction)
     */
    ResultFunction(std::vector<ElementType>* fes, const ResultRequirements& req, UserFunction userFunction = {})
        : gridView{fes->at(0).localView().globalBasis().gridView()},
          resultRequirements_{req},
          fes_{fes},
          userFunction_{userFunction} {
      if constexpr (!std::is_same_v<UserFunction, Impl::DefaultUserFunction>) userFunction_ = userFunction;
    }

  private:
    double evaluateComponent(int eleID, const Dune::FieldVector<ctype, griddim>& local, int comp) const {
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

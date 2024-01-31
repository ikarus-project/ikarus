// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file resultfunction.hh
 * \brief Ikarus Result Function for Stress and other finite element results
 * \ingroup io
 *
 */

#pragma once

#include <type_traits>
#include <utility>

#include <dune/vtk/function.hh>
#include <dune/vtk/vtkwriter.hh>

#include <ikarus/finiteelements/ferequirements.hh>

namespace Ikarus {
namespace Impl {
  struct DefaultUserFunction
  {
  };
} // namespace Impl

/**
 * \brief Wrapper to evaluate results for a vtkwriter.
 * \details
 * Usage:
 * @code
 *   // Usage with Dune::Vtk::VtkWriter
 *   auto resultFunction = Ikarus::ResultFunction<FiniteElement, resType>::asVtkFunction(&fes, feRequirements);
 *   vtkwriter.addPointData(resultFunction);
 *
 *   // Usage with the native Dune::VTKWriter
 *   auto resultFunction = Ikarus::ResultFunction<FiniteElement, resType>::asShared(&fes, feRequirements);
 *   vtkWriter.addVertexData(resultFunction);
 * @endcode
 * @ingroup io
 * @tparam ElementType_ Type of the finite element
 * @tparam resType requested result type
 * @tparam UserFunction Type of the user-defined function for custom result evaluation (default is
DefaultUserFunction)
 */
template <typename FE, ResultType resType, typename UserFunction = Impl::DefaultUserFunction>
class ResultFunction : public Dune::VTKFunction<typename FE::GridView>
{
public:
  using FiniteElement          = FE;
  using FERequirementType      = typename FiniteElement::FERequirementType;
  using GridView               = typename FiniteElement::GridView;
  using ctype                  = typename GridView::ctype;
  constexpr static int griddim = GridView::dimension;
  using Entity                 = typename GridView::template Codim<0>::Entity;

  /**
   * \brief Evaluate the component at a given entity and local coordinates.
   *
   * This function is required by the Dune::VTKFunction interface.
   *
   * \param comp Stress component index
   * \param e Entity on which to evaluate the stress
   * \param local Local coordinates within the entity
   * \return Stress component value
   */
  double evaluate(int comp, const Entity& e, const Dune::FieldVector<ctype, griddim>& local) const override {
    auto index = gridView_.indexSet().index(e);
    return evaluateComponent(index, local, comp);
  }

  /**
   * \brief Get the number of components.
   *
   * This function is required by the Dune::VTKFunction interface.
   *
   * \return Number of stress components
   */
  [[nodiscard]] int ncomps() const override {
    if constexpr (std::is_same_v<UserFunction, Impl::DefaultUserFunction>) {
      Dune::FieldVector<ctype, griddim> val(0.0);

      auto sigma = fes_->at(0).template calculateAt<resType>(feRequirements_, val);

      return static_cast<int>(sigma.rows() * sigma.cols());
    } else
      return userFunction_.ncomps();
  }

  /**
   * \brief Get the name of the result type.
   *
   * This function is required by the Dune::VTKFunction interface.
   *
   * \return String representing the name of the result type
   */
  [[nodiscard]] constexpr std::string name() const override {
    if constexpr (std::is_same_v<UserFunction, Impl::DefaultUserFunction>)
      return toString(resType);
    else
      return userFunction_.name();
  }

  /**
   * \brief Constructor for ResultFunction.
   *
   * Constructs a ResultFunction object with given finite elements, ferequirements
   *
   * \param fes Pointer to a vector of finite elements
   * \param req FERequirements for evaluation
   */
  ResultFunction(std::vector<FiniteElement>* fes, const FERequirementType& req)
      : gridView_{fes->at(0).localView().globalBasis().gridView()},
        feRequirements_{req},
        fes_{fes},
        userFunction_{UserFunction{}} {}

  /**
   * @brief Creates the ResultFunction as a function that can be used with dune-vtk
   *
   * Constructs a ResultFunction object with given finite elements, ferequirements as `Dune::Vtk::Function`
   *
   * @param fes Pointer to a vector of finite elements
   * @param req FERequirements for evaluation
   */
  static auto asVtkFunction(std::vector<ElementType>* fes, const FERequirementType& req) {
    return Dune::Vtk::Function<GridView>(std::make_shared<ResultFunction>(fes, req));
  }

  /**
   * @brief Creates the ResultFunction as a Localfunction that can queried element-wise
   *
   * Constructs a ResultFunction object with given finite elements, ferequirements as `Dune::Vtk::LocalFunction`
   * This leverages the Dune::VTK localfunction interface to be able to query element-wise
   * Usage:
   * @code
   * auto localResultFunction = Ikarus::ResultFunction<FiniteElement, resType>::asLocalFunction(&fes,
   * feRequirements); localResultFunction.bind(element); auto result = localResultFunction(pos);
   * @endcode
   * @param fes Pointer to a vector of finite elements
   * @param req FERequirements for evaluation
   */
  static auto asLocalFunction(std::vector<ElementType>* fes, const FERequirementType& req) {
    return localFunction(asVtkFunction(fes, req));
  }

  /**
   * @brief Creates the ResultFunction as a function that can be used with dune-vtk
   *
   * Constructs a ResultFunction object with given finite elements, ferequirements as shared_ptr to be used with
   * the native Dune VTKWriter
   *
   * @param fes Pointer to a vector of finite elements
   * @param req FERequirements for evaluation
   */
  static auto asShared(std::vector<ElementType>* fes, const FERequirementType& req) {
    return std::make_shared<ResultFunction>(fes, req);
  }

private:
  double evaluateComponent(int eleID, const Dune::FieldVector<ctype, griddim>& local, int comp) const {
    auto result = fes_->at(eleID).template calculateAt<resType>(feRequirements_, local);

    if constexpr (!std::is_same_v<UserFunction, Impl::DefaultUserFunction>)
      return userFunction_(result, comp);
    else
      return result(comp);
  }

  GridView gridView_;
  FERequirementType feRequirements_;
  std::vector<FiniteElement>* fes_;
  [[no_unique_address]] std::string name_{};
  UserFunction userFunction_;
};
} // namespace Ikarus

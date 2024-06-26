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
 * \code
 * // Usage with Dune::Vtk::VtkWriter
 * auto resultFunction = Ikarus::makeResultVtkFunction<resType>(&fes, feRequirements);
 * vtkwriter.addPointData(resultFunction);
 *
 * // Usage with the native Dune::VTKWriter
 * auto resultFunction = Ikarus::makeResultFunction<resType>(&fes, feRequirements);
 * vtkWriter.addVertexData(resultFunction);
 * \endcode
 * \ingroup io
 * \relates
 * \tparam FE Type of the finite element
 * \tparam RT requested result type
 * \tparam UserFunction Type of the user-defined function for custom result evaluation (default is
DefaultUserFunction)
 */
template <typename FE, template <typename, int, int> class RT, typename UserFunction = Impl::DefaultUserFunction>
class ResultFunction : public Dune::VTKFunction<typename FE::GridView>
{
public:
  using FiniteElement          = FE;
  using FERequirementType      = typename FiniteElement::Requirement;
  using GridView               = typename FiniteElement::GridView;
  using ctype                  = typename GridView::ctype;
  constexpr static int griddim = GridView::dimension;
  using Entity                 = typename GridView::template Codim<0>::Entity;

  /**
   * \brief Evaluate the component at a given entity and local coordinates.
   * \details
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
   * \details
   * This function is required by the Dune::VTKFunction interface.
   *
   * \return Number of stress components
   */
  [[nodiscard]] int ncomps() const override {
    if constexpr (std::is_same_v<UserFunction, Impl::DefaultUserFunction>) {
      Dune::FieldVector<ctype, griddim> val(0.0);
      auto sigma = fes_->at(0).template calculateAt<RT>(feRequirements_, val).asVec();
      return static_cast<int>(sigma.size());
    } else
      return userFunction_.ncomps();
  }

  /**
   * \brief Get the name of the result type.
   * \details
   * This function is required by the Dune::VTKFunction interface.
   *
   * \return String representing the name of the result type
   */
  [[nodiscard]] constexpr std::string name() const override {
    if constexpr (requires { userFunction_.name(); })
      return userFunction_.name();
    else
      return toString<RT>();
  }

  /**
   * \brief Get the precision used for this result
   * \details
   * This function is part of the Dune::VTKFunction interface.
   * This has no affect when the ResultFunction is used as part of a the dune-vtk module
   *
   * \return Precision (i.e. float64 or float32)
   */
  Dune::VTK::Precision precision() const override { return prec_; }

  /**
   * \brief Constructor for ResultFunction.
   * \details
   * Constructs a ResultFunction object with given finite elements, ferequirements
   *
   * \param fes Pointer to a vector of finite elements
   * \param req FERequirements for evaluation
   */
  ResultFunction(std::vector<FiniteElement>* fes, const FERequirementType& req,
                 Dune::VTK::Precision prec = Dune::VTK::Precision::float64)
      : gridView_{fes->at(0).localView().globalBasis().gridView()},
        feRequirements_{req},
        fes_{fes},
        prec_{prec},
        userFunction_{UserFunction{}} {}

private:
  double evaluateComponent(int eleID, const Dune::FieldVector<ctype, griddim>& local, int comp) const {
    auto result = fes_->at(eleID).template calculateAt<RT>(feRequirements_, local).asVec();

    if constexpr (!std::is_same_v<UserFunction, Impl::DefaultUserFunction>)
      return userFunction_(result, comp);
    else
      return result(comp);
  }

  GridView gridView_;
  FERequirementType feRequirements_;
  std::vector<FiniteElement>* fes_;
  Dune::VTK::Precision prec_;
  [[no_unique_address]] std::string name_{};
  UserFunction userFunction_;
};

/**
 * \brief Function to create a ResultFunction as a shared_ptr
 * \details
 * Constructs a ResultFunction object with given finite elements, ferequirements as shared_ptr to be used with
 * the native Dune VTKWriter
 *
 * \param fes Pointer to a vector of finite elements
 * \param req FERequirements for evaluation
 * \tparam FE Type of the finite element
 * \tparam RT requested result type
 * \tparam UserFunction Type of the user-defined function for custom result evaluation (default is DefaultUserFunction)
 */
template <template <typename, int, int> class RT, typename UserFunction = Impl::DefaultUserFunction, typename FE>
auto makeResultFunction(std::vector<FE>* fes, const typename FE::Requirement& req,
                        Dune::VTK::Precision prec = Dune::VTK::Precision::float64) {
  return std::make_shared<ResultFunction<FE, RT, UserFunction>>(fes, req, prec);
}

/**
 * \brief Function to create a ResultFunction as a gridfunction that can be used with dune-vtk
 * \details
 * Constructs a ResultFunction object with given finite elements, ferequirements as a VTK::Function to be used with
 * dune-vtk It is possible to construct a localFunction from this as follows
 * \code
 * auto localResultFunction = localFunction(vtkResultFunction);
 * localResultFunction.bind(element);
 * \endcode
 * \param fes Pointer to a vector of finite elements
 * \param req FERequirements for evaluation
 * \tparam FE Type of the finite element
 * \tparam RT requested result type
 * \tparam UserFunction Type of the user-defined function for custom result evaluation (default is
 * DefaultUserFunction)
 */
template <template <typename, int, int> class RT, typename UserFunction = Impl::DefaultUserFunction, typename FE>
auto makeResultVtkFunction(std::vector<FE>* fes, const typename FE::Requirement& req) {
  return Dune::Vtk::Function<typename FE::GridView>(std::make_shared<ResultFunction<FE, RT, UserFunction>>(fes, req));
}

} // namespace Ikarus

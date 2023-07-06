// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_COMPOSEDGRIDFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_COMPOSEDGRIDFUNCTION_HH

#include <type_traits>
#include <tuple>

#include <dune/common/referencehelper.hh>
#include <dune/common/typeutilities.hh>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>


namespace Dune {
namespace Functions {



/**
 * \brief Composition of grid functions with another function.
 *
 * \ingroup FunctionImplementations
 *
 * For given inner grid functions `g0, ..., gn` and an
 * outer function `f` this creates a grid function
 * representing `f(g0(x), ..., gn(x))`. The only assumption
 * made, is that the range types of the inner functions
 * can be passed to the outer ones, and that all grid
 * functions are defined on the same EntitySet.
 *
 * Notice that all functions are captured by value.
 * To store references you can pass `std::ref()`.
 *
 * \tparam OF Type of outer function. std::reference_wrapper is supported.
 * \tparam IF Types of inner outer functions. `std::reference_wrapper` is supported.
 */
template<class OF, class... IF>
class ComposedGridFunctionMod
{
  using InnerFunctions = std::tuple<IF...>;
  using InnerLocalFunctions = std::tuple<decltype(localFunction(resolveRef(std::declval<const IF&>())))...>;

  template<std::size_t i>
  using InnerFunction = std::decay_t<ResolveRef_t<std::tuple_element_t<i, InnerFunctions>>>;

  using OuterFunction = OF;

 public:

  using EntitySet = typename InnerFunction<0>::EntitySet;
  using Element = typename EntitySet::Element;

  using Domain = typename EntitySet::GlobalCoordinate;
  using LocalDomain = Dune::FieldVector<double,3>;

  using Range = decltype(std::declval<OF>()(std::declval<IF>()(std::declval<LocalDomain>())...,double()));

 private:

  using Traits = Imp::GridFunctionTraits<Range(LocalDomain), EntitySet, DefaultDerivativeTraits, 16>;

  class LocalFunction
  {
   public:
    /**
     * \brief Construct the local-function.
     *
     * The local-functions is created from the outer-function of the
     * grid-function and the local-functions of the stored inner-functions.
     **/
    LocalFunction(const ComposedGridFunctionMod& globalFunction) :
        globalFunction_(globalFunction),
        innerLocalFunctions_(globalFunction.innerLocalFunctions())
    {}

    /**
     * \brief Bind the inner local-functions to an `element`.
     *
     * \b Expects:
     * - The `element` is in the entitySet of all inner local-functions.
     *
     * \b Ensures:
     * - All inner local-functions are bound to the same `element`.
     **/
    void bind(const Element& element)
    {
      std::apply([&](auto&... innerFunction) {
        (innerFunction.bind(element),...);
      }, innerLocalFunctions_);
    }

    //! \brief Unbind the inner local-functions.
    void unbind()
    {
      std::apply([&](auto&... innerFunction) {
        (innerFunction.unbind(),...);
      }, innerLocalFunctions_);
    }

    /** \brief Return if the local function is bound to a grid element
     */
    bool bound() const
    {
      return std::apply([](const auto&... innerFunction) {
        return (innerFunction.bound() && ...);
      }, innerLocalFunctions_);
    }

    /**
     * \brief Evaluation of the composed local-function.
     *
     * Returns the outer-function evaluated with the evaluated inner
     * local-functions. The functions are evaluated in the local coordinates `x`.
     *
     * \b Expects:
     * - All inner local-functions are bound to the same element.
     **/
    Range operator()(const LocalDomain& x) const
    {
      Dune::FieldVector<double,2> coords2D= {x[0],x[1]};
      return std::apply([&](const auto&... innerFunction) {
        return globalFunction_.outerFunction_(innerFunction(coords2D)...,x[2]);
      }, innerLocalFunctions_);
    }

    /**
     * \brief Return the local context all inner local-functions are bound to.
     *
     * Since all inner local-functions are bound to the same element, return
     * the local context of either of the inner local-functions.
     *
     * \b Requirements:
     * - Number of inner local-functions > 0
     **/
    const Element& localContext() const
    {
      return std::get<0>(innerLocalFunctions_).localContext();
    }

    //! Not implemented
    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

   private:
    const ComposedGridFunctionMod& globalFunction_;
    InnerLocalFunctions innerLocalFunctions_;
  };

 public:

  /**
   * \brief Create ComposedGridFunctionMod.
   *
   * Outer and inner functions will be captured by value.
   * To store references you can pass `std::ref()`.
   *
   * \param outerFunction The outer function to be composed with the grid functions.
   * \param innerFunctions The inner grid functions
   */
  template<class OFT, class... IFT,
      disableCopyMove<ComposedGridFunctionMod, OFT> = 0,
      std::enable_if_t<(sizeof...(IFT) > 0), int> = 0>
  ComposedGridFunctionMod(OFT&& outerFunction, IFT&&... innerFunctions) :
      outerFunction_(std::forward<OFT>(outerFunction)),
      innerFunctions_(std::forward<IFT>(innerFunctions)...)
  {}

  //! Evaluation of the composed grid function in coordinates `x`
  Range operator()(const Domain& x) const
  {
    return std::apply([&](const auto&... innerFunction) {
      return outerFunction_(innerFunction(x)...);
    }, innerFunctions_);
  }

  //! Not implemented.
  friend typename Traits::DerivativeInterface derivative(const ComposedGridFunctionMod& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  /**
   * \brief Create a local-function of this composed grid-function.
   *
   * The LocalFunction is defined by composition of the outer-function with
   * the corresponding local-functions of the inner-functions.
   **/
  friend LocalFunction localFunctionMod(const ComposedGridFunctionMod& cgf)
  {
    return LocalFunction(cgf);
  }

  /**
   * \brief Return the EntitySet associated to this composed grid-function.
   *
   * It is implicitly assumed that all inner-functions can be bound to the
   * same set of entities. Thus, the EntitySet is either of the EntitySets
   * of the inner-function, e.g., the one from the first inner-function.
   *
   * \b Requirements:
   * - Number of inner-functions > 0
   **/
  const EntitySet& entitySet() const
  {
    return resolveRef(std::get<0>(innerFunctions_)).entitySet();
  }

 protected:

  InnerLocalFunctions innerLocalFunctions() const
  {
    return std::apply([&](const auto&... innerFunction) {
      return std::make_tuple(localFunction(resolveRef(innerFunction))...);
    }, innerFunctions_);
  }

  OuterFunction outerFunction_;
  InnerFunctions innerFunctions_;
};

// deduction guides
template<class OF, class... IF>
ComposedGridFunctionMod(const OF&, const IF&...)
-> ComposedGridFunctionMod<OF,IF...>;


/**
 * \brief Create a ComposedGridFunctionMod that composes grid-functions with another function.
 *
 * \ingroup FunctionImplementations
 *
 * For given inner grid-functions `g0, ..., gn` and an
 * outer-function `f` this creates a grid-function
 * representing `f(g0(x), ..., gn(x))`. The only assumption
 * made, is that the range types of the inner-functions
 * can be passed to the outer ones, and that all grid-functions
 * are defined on the same `EntitySet`.
 *
 * Notice that all functions are captured by value.
 * To store references you can pass `std::ref()`.
 *
 * \param outerFunction The outer-function to be composed with the grid-functions.
 * \param innerFunctions The inner grid-functions
 *
 * \returns A grid-function defined on the same `EntitySet` as the input-functions.
 *
 * \relatesalso ComposedGridFunctionMod
 */
template<class OF, class... IF>
auto makeComposedGridFunctionMod(OF&& outerFunction, IF&&... innerFunction)
{
  using ComposedGridFunctionModType = ComposedGridFunctionMod<std::decay_t<OF>, std::decay_t<IF>...>;
  return ComposedGridFunctionModType(std::forward<OF>(outerFunction), std::forward<IF>(innerFunction)...);
}



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_COMPOSEDGRIDFUNCTION_HH

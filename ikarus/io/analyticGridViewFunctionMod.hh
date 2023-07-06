// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH

#include <type_traits>
#include <optional>

#include <dune/common/typeutilities.hh>

#include <dune/functions/common/signature.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/localderivativetraits.hh>


namespace Dune {
namespace Functions {

namespace Imp {

template<class Signature, class GV, class FLocal, template<class> class DerivativeTraits=DefaultDerivativeTraits>
class LocalAnalyticGridViewFunctionMod;
template <typename VectorType>
VectorType cross(const VectorType& a, const VectorType& b) {

  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
template<class Range, class LocalDomain, class GV, class F, template<class> class DerivativeTraits>
class LocalAnalyticGridViewFunctionMod<Range(LocalDomain), GV, F, DerivativeTraits>
{
 public:
  using Signature = Range(LocalDomain);
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(LocalDomain);

  using GridView = GV;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;
//  using Geometry = typename Element::Geometry;
  using Geometry = typename std::decay<typename Element::Geometry>::type;

  // Use the indirection via derivativeIfImplemented to also support
  // function types F that do not implement derivative. In this case
  // the interface type DifferentiableFunction is using a dummy for
  // the derivative type
  using DerivativeDummy = DifferentiableFunction<DerivativeSignature>;
  using GlobalRawDerivative = decltype(Imp::derivativeIfImplemented<DerivativeDummy, F>(std::declval<F>()));
  using LocalDerivative = LocalAnalyticGridViewFunctionMod<DerivativeSignature, GridView, GlobalRawDerivative, DerivativeTraits>;
  double thickness;
  //! Create the local-function by storing the mapping `f` by value
  template<class FT, disableCopyMove<LocalAnalyticGridViewFunctionMod, FT> = 0>
  LocalAnalyticGridViewFunctionMod(FT&& f,double p_thickness) :
      f_(std::forward<FT>(f))
      ,thickness{p_thickness}
  {}

  //! Constructor that copies the state of the passed element and geometry
  template<class FT>
  LocalAnalyticGridViewFunctionMod(FT&& f, const Element& element, const std::optional<Geometry>& geometry) :
      f_(std::forward<FT>(f)),
      element_(element),
      geometry_(geometry)
  {}


  /**
   * \brief Bind the local-function to an `element`.
   *
   * Stores a copy of the `element` and of its geometry in the class.
   *
   * \b Ensures:
   * - Local-function is bound to the `element`.
   **/
  void bind(const Element& element)
  {
    element_ = element;
    geometry_.emplace(element_.geometry());
  }

  //! Release the bound geometry
  void unbind()
  {
    geometry_.reset();
  }

  /** \brief Return if the local function is bound to a grid element
   */
  bool bound() const
  {
    return static_cast<bool>(geometry_);
  }

  /**
   * \brief Evaluate the stored function `f` in global coordinates mapped by the geometry.
   *
   * \b Expects:
   * - The local-function is bound to an element in bind().
   *
   * \param x  Local coordinate in the bound element
   * \return Evaluation of `f` in global coordinates `x`
   **/
  Range operator()(const LocalDomain& x) const
  {
    assert(!!geometry_);
    Dune::FieldVector<double,2> coords2D= {x[0],x[1]};
    auto globalCoords= geometry_->global(coords2D);
    auto J= geometry_->jacobianTransposed(coords2D);
    auto globalnormal= cross(J[0],J[1]);
    globalnormal/=globalnormal.two_norm();
    return f_(globalCoords+(2*x[2]-1)*globalnormal*thickness/2.0);
  }

  //! Return the bound element
  const Element& localContext() const
  {
    assert(!!geometry_);
    return element_;
  }

  /**
   * \brief Return a local-function representing the derivative.
   *
   * This function computes the derivative of the wrapped function `f`, if
   * available, otherwise use a dummy representation. If the local-function
   * was bound to an element so is its derivative. Otherwise it must be bound
   * before it can be evaluated.
   **/
  friend LocalDerivative derivative(const LocalAnalyticGridViewFunctionMod& t)
  {
    return LocalDerivative(Imp::derivativeIfImplemented<DerivativeDummy, F>(t.f_), t.element_, t.geometry_);
  }

 private:
  F f_;
  Element element_;
  std::optional<Geometry> geometry_ = std::nullopt;
};

} // end namespace Imp




template<class Signature, class GV, class F, template<class> class DerivativeTraits=DefaultDerivativeTraits>
class AnalyticGridViewFunctionMod;


/**
 * \brief Class wrapping any differentiable function as grid function.
 *
 * \ingroup FunctionImplementations
 */
template<class Range, class Domain, class GV, class F, template<class> class DerivativeTraits>
class AnalyticGridViewFunctionMod<Range(Domain), GV, F, DerivativeTraits>
{
 public:
  using Signature = Range(Domain);
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

  using GridView = GV;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;
  using Geometry = typename Element::Geometry;

  // Use the indirection via derivativeIfImplemented to also support
  // function types F that do not implement derivative. In this case
  // the interface type DifferentiableFunction is used a dummy for
  // the derivative type
  using DerivativeDummy = DifferentiableFunction<DerivativeSignature>;
  using GlobalRawDerivative = decltype(Imp::derivativeIfImplemented<DerivativeDummy, F>(std::declval<F>()));
  using Derivative = AnalyticGridViewFunctionMod<DerivativeSignature, GridView, GlobalRawDerivative, DerivativeTraits>;

  using LocalDomain = Dune::FieldVector<double,3>;
  using LocalFunction = typename Imp::LocalAnalyticGridViewFunctionMod<Range(LocalDomain), GridView, F, LocalDerivativeTraits<EntitySet, DerivativeTraits>::template Traits>;
  double thickness;
  //! Create the grid-function by wrapping a function `f` and create a GridViewEntitySet.
  template<class FT>
  AnalyticGridViewFunctionMod(FT&& f, const GridView& gridView,double p_thickness) :
      f_(std::forward<FT>(f)),
      entitySet_(gridView),
      thickness{p_thickness}
  {}

  //! Evaluate the wrapped function `f` directly in global coordinates `x`.
  Range operator()(const Domain& x) const
  {
    return f_(x);
  }

  //! Create a derivative grid-function by wrapping the derivative of `f`.
  friend Derivative derivative(const AnalyticGridViewFunctionMod& t)
  {
    return Derivative(Imp::derivativeIfImplemented<DerivativeDummy, F>(t.f_), t.entitySet_.gridView());
  }

  //! Construct the associated local-function.
  friend LocalFunction localFunctionMod(const AnalyticGridViewFunctionMod& t)
  {
    return LocalFunction(t.f_,t.thickness);
  }

  //! Return the set of entities this local-function can be bound to.
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

 private:
  F f_;
  EntitySet entitySet_;
};



/**
 * \brief Create an AnalyticGridViewFunctionMod from a function and a grid view.
 *
 * \ingroup FunctionImplementations
 *
 * The returned function supports `localFunction()` and stores a copy of the
 * original function.
 *
 * \param f A function object supporting evaluation with global coordinates
 *          of the passed `gridView`.
 * \param gridView The GridView the function should act on.
 *
 * \returns A function that models the GridFunction interface.
 *
 * \relatesalso AnalyticGridViewFunctionMod
 */
template<class F, class GridView>
AnalyticGridViewFunctionMod<
    typename std::invoke_result<F, typename GridView::template Codim<0>::Geometry::GlobalCoordinate>::type  // Range
        (typename GridView::template Codim<0>::Geometry::GlobalCoordinate),                                 // Domain
    GridView,
    typename std::decay<F>::type >                                                                      // Raw type of F (without & or &&)
makeAnalyticGridViewFunctionMod(F&& f, const GridView& gridView,double thickness)
{
  using Domain = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
  using Range = typename std::invoke_result<F, Domain>::type;
  using FRaw = typename std::decay<F>::type;

  return AnalyticGridViewFunctionMod<Range(Domain), GridView, FRaw>(std::forward<F>(f), gridView,thickness);
}



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_ANALYTICGRIDVIEWFUNCTION_HH

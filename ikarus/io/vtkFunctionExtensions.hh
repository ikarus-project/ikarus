
#pragma once

#include <numeric>
#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/common/version.hh>

#include <dune/vtk/localfunction.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/arguments.hh>
#include <dune/vtk/utility/concepts.hh>
#include <dune/vtk/vtkwriterinterface.hh>


#include <memory>
#include <type_traits>

#include <dune/common/typetraits.hh>

//#include "localfunctioninterface.hh"
//#include "legacyvtkfunction.hh"
//#include "defaultvtkfunction.hh"

#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>



#include "vtkfunctionmod.hh"


namespace Dune
{
namespace Vtk
{
/// \brief An abstract base class for LocalFunctionMods that can be bound to an element and
/// evaluated in local coordinates w.r.t. to a component of its value.
template <class GridView>
class LocalFunctionModInterface
{
 public:
  using Entity = typename GridView::template Codim<0>::Entity;
  using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

  /// Bind the function to the grid entity
  virtual void bind (Entity const& entity) = 0;

  /// Unbind from the currently bound entity
  virtual void unbind () = 0;

  /// Evaluate single component comp in the entity at local coordinates xi
//  virtual double evaluate (int comp, LocalCoordinate const& xi) const = 0;
  virtual double evaluate (int comp, Dune::FieldVector<double,3> const& xi) const = 0;

  /// Virtual destructor
  virtual ~LocalFunctionModInterface () = default;
};

} // end namespace Vtk
} // end namespace DuneFunctionM
namespace Dune
{
namespace Vtk
{
/// Type erasure for dune-functions LocalFunction interface
template <class GridView, class LocalFunction>
class LocalFunctionWrapperMod final
    : public LocalFunctionModInterface<GridView>
{
  using Self = LocalFunctionWrapperMod;
  using Interface = LocalFunctionModInterface<GridView>;
  using Entity = typename Interface::Entity;
  using LocalCoordinate = typename Interface::LocalCoordinate;

 public:
  /// Constructor. Stores a copy of the passed `localFct` in a local variable.
  template <class LocalFct,
      disableCopyMove<Self, LocalFct> = 0>
  explicit LocalFunctionWrapperMod (LocalFct&& localFct)
      : localFct_(std::forward<LocalFct>(localFct))
  {}

  /// Bind the LocalFunction to the Entity
  virtual void bind (Entity const& entity) override
  {
    localFct_.bind(entity);
  }

  /// Unbind the LocalFunction from the Entity
  virtual void unbind () override
  {
    localFct_.unbind();
  }

  /// Evaluate the LocalFunction in LocalCoordinates
//  virtual double evaluate (int comp, LocalCoordinate const& xi) const override
//  {
//    return evaluateImpl(comp, localFct_(xi));
//  }

  virtual double evaluate (int comp, Dune::FieldVector<double,3> const& xi) const override
  {
    return evaluateImpl(comp, localFct_(xi));
  }

 private:
  // Evaluate a component of a vector valued data
  template <class T, int N, int M>
  double evaluateImpl (int comp, FieldMatrix<T,N,M> const& mat) const
  {
    int r = comp / 3;
    int c = comp % 3;
    return r < N && c < M ? mat[r][c] : 0.0;
  }

  // Evaluate a component of a vector valued data
  template <class T, int N>
  double evaluateImpl (int comp, FieldVector<T,N> const& vec) const
  {
    return comp < N ? vec[comp] : 0.0;
  }

  // Evaluate a component of a vector valued data
  template <class T,
      std::enable_if_t<IsIndexable<T,int>::value, int> = 0>
  double evaluateImpl (int comp, T const& value) const
  {
    return value[comp];
  }

  // Return the scalar values
  template <class T,
      std::enable_if_t<not IsIndexable<T,int>::value, int> = 0>
  double evaluateImpl (int comp, T const& value) const
  {
    assert(comp == 0);
    return value;
  }

 private:
  LocalFunction localFct_;
};

} // end namespace Vtk
} // end namespace Dune

#include <memory>

#include <dune/grid/io/file/vtk/function.hh>


namespace Dune
{
namespace Vtk
{
/// Type erasure for Legacy VTKFunction
template <class GridView>
class VTKLocalFunctionModWrapper final
    : public LocalFunctionModInterface<GridView>
{
  using Interface = LocalFunctionModInterface<GridView>;
  using Entity = typename Interface::Entity;
  using LocalCoordinate = typename Interface::LocalCoordinate;

 public:
  /// Constructor. Stores a shared pointer to the passed Dune::VTKFunction
  explicit VTKLocalFunctionModWrapper (std::shared_ptr<Dune::VTKFunctionMod<GridView> const> const& fct)
      : fct_(fct)
  {}

  /// Stores a pointer to the passed entity
  virtual void bind (Entity const& entity) override
  {
    entity_ = &entity;
  }

  /// Unsets the stored entity pointer
  virtual void unbind () override
  {
    entity_ = nullptr;
  }

  /// Evaluate the Dune::VTKFunction in LocalCoordinates on the stored Entity
  virtual double evaluate (int comp, Dune::FieldVector<double,3> const& xi) const override
  {
    return fct_->evaluate(comp, *entity_, xi);
  }

 private:
  std::shared_ptr<VTKFunctionMod<GridView> const> fct_;
  Entity const* entity_;
};

} //end namespace Vtk
} // end namespace Dune


namespace Dune
{
namespace Vtk
{
/// \brief A Vtk::LocalFunctionMod is a function-like object that can be bound to a grid element
/// an that provides an evaluate method with a component argument.
/**
* Stores internally a Vtk::LocalFunctionModInterface object for the concrete evaluation.
**/
template <class GridView>
class LocalFunctionMod
{
  using Self = LocalFunctionMod;
  using Entity = typename GridView::template Codim<0>::Entity;
  using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

  template <class LF, class E>
  using HasBind = decltype((std::declval<LF>().bind(std::declval<E const&>()), true));

 private:
  struct RangeProxy
  {
    using value_type = double;
    using field_type = double;

    RangeProxy (LocalFunctionModInterface<GridView> const& localFct,
                std::vector<int> const& components,
                LocalCoordinate const& local)
        : localFct_(localFct)
        , components_(components)
        , local_(local)
    {}

    std::size_t size () const
    {
      return components_.size();
    }

    double operator[] (std::size_t i) const
    {
      return i < size() ? localFct_.evaluate(components_[i], local_) : 0.0;
    }

   private:
    LocalFunctionModInterface<GridView> const& localFct_;
    std::vector<int> const& components_;
    LocalCoordinate local_;
  };

 public:
  /// Construct the Vtk::LocalFunctionMod from any function object that has a bind(element) method.
  template <class LF,
      disableCopyMove<Self, LF> = 0,
      HasBind<LF,Entity> = true>
  explicit LocalFunctionMod (LF&& lf)
      : localFct_(std::make_shared<LocalFunctionWrapperMod<GridView,LF>>(std::forward<LF>(lf)))
  {}

  /// Construct a Vtk::LocalFunctionMod from a legacy VTKFunction
  explicit LocalFunctionMod (std::shared_ptr<Dune::VTKFunctionMod<GridView> const> const& lf)
      : localFct_(std::make_shared<VTKLocalFunctionModWrapper<GridView>>(lf))
  {}

  /// Allow the default construction of a Vtk::LocalFunctionMod. After construction, the
  /// LocalFunctionMod is in an invalid state.
  LocalFunctionMod () = default;

  /// Bind the function to the grid entity
  void bind (Entity const& entity)
  {
    assert(bool(localFct_));
    localFct_->bind(entity);
  }

  /// Unbind from the currently bound entity
  void unbind ()
  {
    assert(bool(localFct_));
    localFct_->unbind();
  }

  /// Return a proxy object to access the components of the range vector
  RangeProxy operator() (LocalCoordinate const& xi) const
  {
    assert(bool(localFct_));
    return {*localFct_, components_, xi};
  }

  /// Evaluate the `c`th component of the Range value at local coordinate `xi`
//  double evaluate (int c, LocalCoordinate const& xi) const
//  {
//    assert(bool(localFct_));
//    return c < int(components_.size()) ? localFct_->evaluate(components_[c], xi) : 0.0;
//  }

  double evaluate (int c, Dune::FieldVector<double,3> const& xi) const
  {
    assert(bool(localFct_));
    return c < int(components_.size()) ? localFct_->evaluate(components_[c], xi) : 0.0;
  }

  void setComponents (std::vector<int> components)
  {
    components_ = std::move(components);
  }

 private:
  std::shared_ptr<LocalFunctionModInterface<GridView>> localFct_ = nullptr;
  std::vector<int> components_;
};

} // end namespace Vtk
} // end namespace Dune


#pragma once

#include <type_traits>
#include <dune/geometry/type.hh>

namespace Dune
{
namespace Vtk
{
//template <class...> struct CheckTypes {};

//template <class DataCollector, class DC = std::decay_t<DataCollector>>
//using IsDataCollector = decltype((
//    std::declval<DC&>().update(),
//        std::declval<DC>().numPoints(),
//        std::declval<DC>().numCells(),
//        CheckTypes<typename DC::GridView>{},
//        true));
//
//template <class GridView, class GV = std::decay_t<GridView>>
//using IsGridView = decltype((
//    std::declval<GV>().grid(),
//        std::declval<GV>().indexSet(),
//        std::declval<GV>().size(0),
//        std::declval<GV>().size(std::declval<Dune::GeometryType>()),
//        CheckTypes<typename GV::Grid, typename GV::IndexSet>{},
//        true));

template <class GridFunction, class GF = std::decay_t<GridFunction>>
using IsGridFunctionMod = decltype((
    localFunctionMod(std::declval<GF const&>()),
        true));

template <class LocalFunction, class LocalContext,
    class LF = std::decay_t<LocalFunction>,
    class LC = std::decay_t<LocalContext>>
using IsLocalFunctionMod = decltype((
    std::declval<LF&>().bind(std::declval<LC const&>()),
        std::declval<LF&>().unbind(),
        true));

} // end namespace Vtk
} // end namespace Dune


namespace Dune
{
// forward declarations
template <class T, int N>
class FieldVector;

template <class T, int N, int M>
class FieldMatrix;

namespace Vtk
{
/// Wrapper class for functions allowing local evaluations.
template <class GridView>
class FunctionMod
{
  using Element = typename GridView::template Codim<0>::Entity;
  using LocalDomain = typename Element::Geometry::LocalCoordinate;

  template <class F, class D>
  using Range = std::decay_t<std::result_of_t<F(D)>>;

 private:

  template <class T, int N>
  static auto sizeOfImpl (FieldVector<T,N>) -> std::integral_constant<int, N> { return {}; }

  template <class T, int N, int M>
  static auto sizeOfImpl (FieldMatrix<T,N,M>) -> std::integral_constant<int, N*M> { return {}; }

  static auto sizeOfImpl (...) -> std::integral_constant<int, 1> { return {}; }

  template <class T>
  static constexpr int sizeOf () { return decltype(sizeOfImpl(std::declval<T>()))::value; }

  static std::vector<int> allComponents(int n)
  {
    std::vector<int> components(n);
    std::iota(components.begin(), components.end(), 0);
    return components;
  }

 public:
  /// (1) Construct from a LocalFunctionMod directly
  /**
  * \param localFct    A local-function, providing a `bind(Element)` and an `operator()(LocalDomain)`
  * \param name        The name to use as identification in the VTK file
  * \param components  A vector of component indices to extract from the range type
  * \param category    The \ref Vtk::RangeTypes category for the range. [Vtk::RangeTypes::UNSPECIFIED]
  * \param dataType    The \ref Vtk::DataTypes used in the output. [Vtk::DataTypes::FLOAT64]
  *
  * The arguments `category` and `dataType` can be passed in any order.
  *
  * NOTE: Stores the localFunctionMod by value.
  **/
  template <class LF, class... Args,
      IsLocalFunctionMod<LF,Element> = true
  >
  FunctionMod (LF&& localFct, std::string name, std::vector<int> components, Args&&... args)
      : localFct_(std::forward<LF>(localFct))
      , name_(std::move(name))
  {
    std::cout<<"(1)"<<std::endl;
    setComponents(std::move(components));
    setRangeType(getArg<Vtk::RangeTypes>(args..., Vtk::RangeTypes::UNSPECIFIED), components_.size());
    setDataType(getArg<Vtk::DataTypes>(args..., Vtk::DataTypes::FLOAT64));
  }

  /// (2) Construct from a LocalFunctionMod directly
  /**
  * \param localFct   A local-function, providing a `bind(Element)` and an `operator()(LocalDomain)`
  * \param name    The name to use as identification in the VTK file
  * \param ncomps  Number of components of the pointwise data. Is extracted
  *                from the range type of the GridFunctionMod if not given.
  *
  * Forwards all the other parmeters to the constructor (1)
  *
  * NOTE: Stores the localFunctionMod by value.
  **/
  template <class LF, class... Args,
      IsLocalFunctionMod<LF,Element> = true
  >
  FunctionMod (LF&& localFct, std::string name, int ncomps, Args&&... args)
      : FunctionMod(std::forward<LF>(localFct), std::move(name), allComponents(ncomps),
                 std::forward<Args>(args)...)
  {    std::cout<<"(2)"<<std::endl;}

  /// (3) Construct from a LocalFunctionMod directly.
  /**
   * Same as Constructor (1) or (2) but deduces the number of components from
   * the static range type of the local-function. This defaults to 1 of no
   * static size information could be extracted.
   **/
  template <class LF, class... Args,
      class R = Range<LF,LocalDomain> ,
      IsLocalFunctionMod<LF,Element> = true
  >
  FunctionMod (LF&& localFct, std::string name, Args&&... args)
      : FunctionMod(std::forward<LF>(localFct), std::move(name), sizeOf<R>(),
                 std::forward<Args>(args)...)
  {    std::cout<<"(3)"<<std::endl;}

  /// (4) Construct from a Vtk::FunctionMod
  template <class... Args>
  explicit FunctionMod (FunctionMod<GridView> const& fct, Args&&... args)
      : FunctionMod(fct.localFct_,
                 getArg<std::string, char const*>(args..., fct.name_),
                 getArg<int,unsigned int,long,unsigned long,std::vector<int>>(args..., fct.components_),
                 getArg<Vtk::RangeTypes>(args..., fct.rangeType_),
                 getArg<Vtk::DataTypes>(args..., fct.dataType_))
  {    std::cout<<"(4)"<<std::endl;}

  /// (5) Construct from a GridFunctionMod
  /**
  * \param fct   A Grid(View)-function, providing a `localFunctionMod(fct)`
  * \param name  The name to use as identification in the VTK file
  *
  * Forwards all other arguments to the constructor (1) or (2).
  *
  * NOTE: Stores the localFunctionMod(fct) by value.
  */
  template <class GF, class... Args,
      disableCopyMove<FunctionMod, GF> = 0, IsGridFunctionMod<GF> = true>
  FunctionMod (GF&& fct, std::string name, Args&&... args)
      : FunctionMod(localFunctionMod(std::forward<GF>(fct)), std::move(name), std::forward<Args>(args)...)
  {    std::cout<<"(5)"<<std::endl;}

  /// (6) Constructor that forwards the number of components and data type to the other constructor
  template <class F>
  FunctionMod (F&& fct, Vtk::FieldInfo info, ...)
      : FunctionMod(std::forward<F>(fct), info.name(), info.size(), info.rangeType(), info.dataType())
  {    std::cout<<"(6)"<<std::endl;}

  /// (7) Automatically extract name and num components from GridFunctionMod if available
  template <class F, class... Args,
      disableCopyMove<FunctionMod, F> = 0,
      IsGridFunctionMod<F> = true,
      class = decltype(std::declval<F>().name()),
      class = decltype(std::declval<F>().numComponents()),
      class = decltype(std::declval<F>().dataType()) >
  explicit FunctionMod (F&& fct, ...)
      : FunctionMod(localFunctionMod(std::forward<F>(fct)), fct.name(), fct.numComponents(),
                 Vtk::RangeTypes::UNSPECIFIED, fct.dataType())
  {    std::cout<<"(7)"<<std::endl;}

  /// (8) Construct from legacy VTKFunctionMod
  /**
  * \param fct  The Dune::VTKFunctionMod to wrap
  **/
  explicit FunctionMod (std::shared_ptr<Dune::VTKFunctionMod<GridView> const> const& fct, ...)
      : localFct_(fct)
      , name_(fct->name())
  {
    std::cout<<"(8)"<<std::endl;
    setComponents(fct->ncomps());
    setDataType(dataTypeOf(fct->precision()));
    setRangeType(rangeTypeOf(fct->ncomps()));
  }

  /// (9) Default constructor. After construction, the function is an an invalid state.
  FunctionMod () = default;

  /// Create a LocalFunctionMod
  friend Vtk::LocalFunctionMod<GridView> localFunctionMod (FunctionMod const& self)
  {
    return self.localFct_;
  }

  /// Return a name associated with the function
  std::string const& name () const
  {
    return name_;
  }

  /// Set the function name
  void setName (std::string name)
  {
    name_ = std::move(name);
  }

  /// Return the number of components of the Range as it is written to the file
  int numComponents () const
  {
    return rangeType_ == Vtk::RangeTypes::SCALAR ? 1 :
           rangeType_ == Vtk::RangeTypes::VECTOR ? 3 :
           rangeType_ == Vtk::RangeTypes::TENSOR ? 9 : int(components_.size());
  }

  /// Set the components of the Range to visualize
  void setComponents (std::vector<int> components)
  {
    components_ = components;
    localFct_.setComponents(components_);
  }

  /// Set the number of components of the Range and generate component range [0...ncomps)
  void setComponents (int ncomps)
  {
    setComponents(allComponents(ncomps));
  }

  /// Return the VTK Datatype associated with the functions range type
  Vtk::DataTypes dataType () const
  {
    return dataType_;
  }

  /// Set the data-type for the components
  void setDataType (Vtk::DataTypes type)
  {
    dataType_ = type;
  }

  /// The category of the range, SCALAR, VECTOR, TENSOR, or UNSPECIFIED
  Vtk::RangeTypes rangeType () const
  {
    return rangeType_;
  }

  /// Set the category of the range, SCALAR, VECTOR, TENSOR, or UNSPECIFIED
  void setRangeType (Vtk::RangeTypes type, std::size_t ncomp = 1)
  {
    rangeType_ = type;
    if (type == Vtk::RangeTypes::AUTO)
      rangeType_ = rangeTypeOf(ncomp);
  }

  /// Set all the parameters from a FieldInfo object
  void setFieldInfo (Vtk::FieldInfo info)
  {
    setName(info.name());
    setComponents(info.size());
    setRangeType(info.rangeType());
    setDataType(info.dataType());
  }

 private:
  Vtk::LocalFunctionMod<GridView> localFct_;
  std::string name_;
  std::vector<int> components_;
  Vtk::DataTypes dataType_ = Vtk::DataTypes::FLOAT64;
  Vtk::RangeTypes rangeType_ = Vtk::RangeTypes::UNSPECIFIED;
};

} // end namespace Vtk
} // end namespace Dune

#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/vtk/filewriter.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/types.hh>

namespace Dune
{
/// Interface for file writers for the Vtk XML file formats
/**
 * \tparam GV  Model of Dune::GridView
 * \tparam DC  Model of \ref DataCollectorInterface
 **/
template <class GV, class DC>
class VtkWriterInterfaceMod
    : public Vtk::FileWriter
{
  template <class> friend class TimeseriesWriter;
  template <class> friend class PvdWriter;

 public:
  using GridView = GV;
  using DataCollector = DC;

 protected:
  using VtkFunction = Dune::Vtk::FunctionMod<GridView>;
  using pos_type = typename std::ostream::pos_type;

  enum PositionTypes {
    POINT_DATA,
    CELL_DATA
  };

 public:
  /// \brief Constructor, passes the gridView to the DataCollector
  /**
   * Creates a new VtkWriterInterfaceMod for the provided GridView. Initializes a
   * DataCollector that is used to collect point coordinates, cell connectivity and
   * data values.
   *
   * This constructor assumes, that the DataCollector can be constructed from a single argument,
   * the passed gridView.
   *
   * \param gridView  Implementation of Dune::GridView
   * \param format    Format of the VTK file, either Vtk::FormatTypes::BINARY, Vtk::FormatTypes::ASCII, or Vtk::COMPRESSED
   * \param datatype  Data type of a single component of the point coordinates [Vtk::DataTypes::FLOAT32]
   * \param headertype  Integer type used in binary data headers [Vtk::DataTypes::UINT32]
   **/
  VtkWriterInterfaceMod (GridView const& gridView,
                      Vtk::FormatTypes format = Vtk::FormatTypes::BINARY,
                      Vtk::DataTypes datatype = Vtk::DataTypes::FLOAT32,
                      Vtk::DataTypes headertype = Vtk::DataTypes::UINT32)
      : VtkWriterInterfaceMod(std::make_shared<DataCollector>(gridView), format, datatype, headertype)
  {}

  /// \brief Constructor, wraps the passed DataCollector in a non-destroying shared_ptr
  VtkWriterInterfaceMod (DataCollector& dataCollector,
                      Vtk::FormatTypes format = Vtk::FormatTypes::BINARY,
                      Vtk::DataTypes datatype = Vtk::DataTypes::FLOAT32,
                      Vtk::DataTypes headertype = Vtk::DataTypes::UINT32)
      : VtkWriterInterfaceMod(stackobject_to_shared_ptr(dataCollector), format, datatype, headertype)
  {}

  /// \brief Constructor, stores the passed DataCollector
  VtkWriterInterfaceMod (std::shared_ptr<DataCollector> dataCollector,
                      Vtk::FormatTypes format = Vtk::FormatTypes::BINARY,
                      Vtk::DataTypes datatype = Vtk::DataTypes::FLOAT32,
                      Vtk::DataTypes headertype = Vtk::DataTypes::UINT32)
      : dataCollector_(std::move(dataCollector))
  {
    setFormat(format);
    setDatatype(datatype);
    setHeadertype(headertype);
  }


  /// \brief Write the attached data to the file
  /**
   * \param fn   Filename of the VTK file. May contain a directory and any file extension.
   * \param dir  The optional parameter specifies the directory of the partition files for parallel writes.
   *
   * \returns File name that is actually written.
   **/
  virtual std::string write (std::string const& fn, std::optional<std::string> dir = {}) const override;

  /// \brief Attach point data to the writer
  /**
   * Attach a global function to the writer that will be evaluated at grid points
   * (vertices and higher order points). The global function must be
   * assignable to the function wrapper \ref Vtk::Function. Additional argument
   * for output datatype and number of components can be passed. See \ref Vtk::Function
   * Constructor for possible arguments.
   *
   * \param fct     A GridFunction, LocalFunction, or Dune::VTKFunction
   * \param args... Additional arguments, like `name`, `numComponents`, `dataType` or `Vtk::FieldInfo`
   **/
  template <class Function, class... Args>
  VtkWriterInterfaceMod& addPointData (Function&& fct, Args&&... args)
  {
    pointData_.emplace_back(std::forward<Function>(fct), std::forward<Args>(args)...,
                            datatype_, Vtk::RangeTypes::AUTO);
    return *this;
  }

  /// \brief Attach cell data to the writer
  /**
   * Attach a global function to the writer that will be evaluated at cell centers.
   * The global function must be assignable to the function wrapper \ref Vtk::Function.
   * Additional argument for output datatype and number of components can be passed.
   * See \ref Vtk::Function Constructor for possible arguments.
   *
   * \param fct     A GridFunction, LocalFunction, or Dune::VTKFunction
   * \param args... Additional arguments, like `name`, `numComponents`, `dataType` or `Vtk::FieldInfo`
   **/
  template <class Function, class... Args>
  VtkWriterInterfaceMod& addCellData (Function&& fct, Args&&... args)
  {
    cellData_.emplace_back(std::forward<Function>(fct), std::forward<Args>(args)...,
                           datatype_, Vtk::RangeTypes::AUTO);
    return *this;
  }


  // Sets the VTK file format
  void setFormat (Vtk::FormatTypes format)
  {
    format_ = format;

    if (format_ == Vtk::FormatTypes::COMPRESSED) {
#if HAVE_VTK_ZLIB
      compressor_ = Vtk::CompressorTypes::ZLIB;
#else
      std::cout << "Dune is compiled without compression. Falling back to BINARY VTK output!\n";
      format_ = Vtk::FormatTypes::BINARY;
#endif
    } else {
      compressor_ = Vtk::CompressorTypes::NONE;
    }

  }

  /// Sets the global datatype used for coordinates and other global float values
  void setDatatype (Vtk::DataTypes datatype)
  {
    datatype_ = datatype;
  }

  /// Sets the integer type used in binary data headers
  void setHeadertype (Vtk::DataTypes datatype)
  {
    headertype_ = datatype;
  }

  /// Sets the compressor type used in binary data headers, Additionally a compression
  /// level can be passed with level = -1 means: default compression level. Level must be in [0-9]
  void setCompressor (Vtk::CompressorTypes compressor, int level = -1)
  {
    compressor_ = compressor;
    compression_level = level;
    VTK_ASSERT(level >= -1 && level <= 9);

    if (compressor_ != Vtk::CompressorTypes::NONE)
      format_ = Vtk::FormatTypes::COMPRESSED;
  }

 private:
  /// Write a serial VTK file in Unstructured format
  virtual void writeSerialFile (std::ofstream& out) const = 0;

  /// Write a parallel VTK file `pfilename.pvtx` in XML format,
  /// with `size` the number of pieces and serial files given by `pfilename_p[i].vtu`
  /// for [i] in [0,...,size).
  virtual void writeParallelFile (std::ofstream& out, std::string const& pfilename, int size) const = 0;

  /// Return the file extension of the serial file (not including the dot)
  virtual std::string fileExtension () const = 0;

  /// Write points and cells in raw/compressed format to output stream
  virtual void writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const = 0;

 protected:
  // Write the point or cell values given by the grid function `fct` to the
  // output stream `out`. In case of binary format, append the streampos of XML
  // attributes "offset" to the vector `offsets`.
  void writeData (std::ofstream& out,
                  std::vector<pos_type>& offsets,
                  VtkFunction const& fct,
                  PositionTypes type,
                  std::optional<std::size_t> timestep = {}) const;

  // Write point-data and cell-data in raw/compressed format to output stream
  void writeDataAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const;

  // Write the coordinates of the vertices to the output stream `out`. In case
  // of binary format, appends the streampos of XML attributes "offset" to the
  // vector `offsets`.
  void writePoints (std::ofstream& out,
                    std::vector<pos_type>& offsets,
                    std::optional<std::size_t> timestep = {}) const;

  // Write Appended section and fillin offset values to XML attributes
  void writeAppended (std::ofstream& out, std::vector<pos_type> const& offsets) const;

  // Write the `values` in blocks (possibly compressed) to the output
  // stream `out`. Return the written block size.
  template <class HeaderType, class FloatType>
  std::uint64_t writeValuesAppended (std::ofstream& out, std::vector<FloatType> const& values) const;

  // Write the `values` in a space and newline separated list of ascii representations.
  // The precision is controlled by the datatype and numerical_limits::digits10.
  template <class T>
  void writeValuesAscii (std::ofstream& out, std::vector<T> const& values) const;

  // Write the XML file header of a VTK file `<VTKFile ...>`
  void writeHeader (std::ofstream& out, std::string const& type) const;

  /// Return PointData/CellData attributes for the name of the first scalar/vector/tensor DataArray
  std::string getNames (std::vector<VtkFunction> const& data) const;

  // Returns endianness
  std::string getEndian () const
  {
    short i = 1;
    return (reinterpret_cast<char*>(&i)[1] == 1 ? "BigEndian" : "LittleEndian");
  }

  // provide accessor to \ref fileExtension virtual method
  std::string getFileExtension () const
  {
    return fileExtension();
  }

  // Returns the VTK file format initialized in the constructor
  Vtk::FormatTypes getFormat () const
  {
    return format_;
  }

  // Returns the global datatype used for coordinates and other global float values
  Vtk::DataTypes getDatatype () const
  {
    return datatype_;
  }

  // Return the global MPI communicator.
  auto comm () const
  {
    return MPIHelper::getCommunication();
  }

 protected:
  std::shared_ptr<DataCollector> dataCollector_;

  Vtk::FormatTypes format_;
  Vtk::DataTypes datatype_;
  Vtk::DataTypes headertype_;
  Vtk::CompressorTypes compressor_ = Vtk::CompressorTypes::NONE;

  // attached data
  std::vector<VtkFunction> pointData_;
  std::vector<VtkFunction> cellData_;

  std::size_t const block_size = 1024*32;
  int compression_level = -1; // in [0,9], -1 ... use default value
};


template <class Writer>
struct IsVtkWriterMod
{
  template <class GV, class DC>
  static std::uint16_t test(VtkWriterInterfaceMod<GV,DC> const&);
  static std::uint8_t  test(...); // fall-back overload

  static constexpr bool value = sizeof(test(std::declval<Writer>())) > sizeof(std::uint8_t);
};

} // end namespace Dune


#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#if HAVE_VTK_ZLIB
#include <zlib.h>
#endif

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/errors.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class GV, class DC>
std::string VtkWriterInterfaceMod<GV,DC>
::write (std::string const& fn, std::optional<std::string> dir) const
{
  dataCollector_->update();

  auto p = Vtk::Path(fn);
  auto name = p.stem();
  p.removeFilename();

  Vtk::Path fn_dir = p;
  Vtk::Path data_dir = dir ? Vtk::Path(*dir) : fn_dir;
  Vtk::Path rel_dir = Vtk::relative(data_dir, fn_dir);

  std::string serial_fn = data_dir.string() + '/' + name.string();
  std::string parallel_fn = fn_dir.string() + '/' + name.string();
  std::string rel_fn = rel_dir.string() + '/' + name.string();

  if (comm().size() > 1)
    serial_fn += "_p" + std::to_string(comm().rank());

  std::string outputFilename;

  { // write serial file
    outputFilename = serial_fn + "." + fileExtension();
    std::ofstream serial_out(outputFilename, std::ios_base::ate | std::ios::binary);
    assert(serial_out.is_open());

    serial_out.imbue(std::locale::classic());
    serial_out << std::setprecision(datatype_ == Vtk::DataTypes::FLOAT32
                                    ? std::numeric_limits<float>::max_digits10
                                    : std::numeric_limits<double>::max_digits10);

    writeSerialFile(serial_out);
  }

  if (comm().size() > 1 && comm().rank() == 0) {
    // write parallel filee
    outputFilename = parallel_fn + ".p" + fileExtension();
    std::ofstream parallel_out(outputFilename, std::ios_base::ate | std::ios::binary);
    assert(parallel_out.is_open());

    parallel_out.imbue(std::locale::classic());
    parallel_out << std::setprecision(datatype_ == Vtk::DataTypes::FLOAT32
                                      ? std::numeric_limits<float>::max_digits10
                                      : std::numeric_limits<double>::max_digits10);

    writeParallelFile(parallel_out, rel_fn, comm().size());
  }

  return outputFilename;
}


template <class GV, class DC>
void VtkWriterInterfaceMod<GV,DC>
::writeData (std::ofstream& out, std::vector<pos_type>& offsets,
             VtkFunction const& fct, PositionTypes type,
             std::optional<std::size_t> timestep) const
{
  out << "<DataArray Name=\"" << fct.name() << "\" type=\"" << to_string(fct.dataType()) << "\""
      << " NumberOfComponents=\"" << fct.numComponents() << "\" format=\"" << (format_ == Vtk::FormatTypes::ASCII ? "ascii\"" : "appended\"");
  if (timestep)
    out << " TimeStep=\"" << *timestep << "\"";

  if (format_ == Vtk::FormatTypes::ASCII) {
    out << ">\n";
    if (type == POINT_DATA)
      writeValuesAscii(out, dataCollector_->template pointData<double>(fct));
    else
      writeValuesAscii(out, dataCollector_->template cellData<double>(fct));
    out << "</DataArray>\n";
  } else {
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
void VtkWriterInterfaceMod<GV,DC>
::writePoints (std::ofstream& out, std::vector<pos_type>& offsets,
               std::optional<std::size_t> timestep) const
{
  out << "<DataArray type=\"" << to_string(datatype_) << "\""
      << " NumberOfComponents=\"3\" format=\"" << (format_ == Vtk::FormatTypes::ASCII ? "ascii\"" : "appended\"");
  if (timestep)
    out << " TimeStep=\"" << *timestep << "\"";

  if (format_ == Vtk::FormatTypes::ASCII) {
    out << ">\n";
    writeValuesAscii(out, dataCollector_->template points<double>());
    out << "</DataArray>\n";
  } else {
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
void VtkWriterInterfaceMod<GV,DC>
::writeDataAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const
{
  for (auto const& v : pointData_) {
    Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(v.dataType(), headertype_,
                                                                [&](auto f, auto h) {
                                                                  using F = typename decltype(f)::type;
                                                                  using H = typename decltype(h)::type;
                                                                  blocks.push_back(this->template writeValuesAppended<H>(out, dataCollector_->template pointData<F>(v)));
                                                                });
  }

  for (auto const& v : cellData_) {
    Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(v.dataType(), headertype_,
                                                                [&](auto f, auto h) {
                                                                  using F = typename decltype(f)::type;
                                                                  using H = typename decltype(h)::type;
                                                                  blocks.push_back(this->template writeValuesAppended<H>(out, dataCollector_->template cellData<F>(v)));
                                                                });
  }
}


template <class GV, class DC>
void VtkWriterInterfaceMod<GV,DC>
::writeAppended (std::ofstream& out, std::vector<pos_type> const& offsets) const
{
  if (is_a(format_, Vtk::FormatTypes::APPENDED)) {
    out << "<AppendedData encoding=\"raw\">\n_";
    std::vector<std::uint64_t> blocks;
    writeGridAppended(out, blocks);
    writeDataAppended(out, blocks);
    out << "</AppendedData>\n";
    pos_type appended_pos = out.tellp();

    pos_type offset = 0;
    for (std::size_t i = 0; i < offsets.size(); ++i) {
      out.seekp(offsets[i]);
      out << '"' << offset << '"';
      offset += pos_type(blocks[i]);
    }

    out.seekp(appended_pos);
  }
}


//namespace Impl {
//
//template <class T, std::enable_if_t<(sizeof(T)>1), int> = 0>
//inline T const& printable (T const& t) { return t; }
//
//inline std::int16_t printable (std::int8_t c) { return std::int16_t(c); }
//inline std::uint16_t printable (std::uint8_t c) { return std::uint16_t(c); }
//
//} // end namespace Impl


template <class GV, class DC>
template <class T>
void VtkWriterInterfaceMod<GV,DC>
::writeValuesAscii (std::ofstream& out, std::vector<T> const& values) const
{
  assert(is_a(format_, Vtk::FormatTypes::ASCII) && "Function should by called only in ascii mode!\n");
  std::size_t i = 0;
  for (auto const& v : values)
    out << Impl::printable(v) << (++i % 6 != 0 ? ' ' : '\n');
  if (i % 6 != 0)
    out << '\n';
}

template <class GV, class DC>
void VtkWriterInterfaceMod<GV,DC>
::writeHeader (std::ofstream& out, std::string const& type) const
{
  out << "<VTKFile"
      << " type=\"" << type << "\""
      << " version=\"1.0\""
      << " header_type=\"" << to_string(headertype_) << "\"";
  if (format_ != Vtk::FormatTypes::ASCII)
    out << " byte_order=\"" << getEndian() << "\"";
  if (compressor_ != Vtk::CompressorTypes::NONE)
    out << " compressor=\"" << to_string(compressor_) << "\"";
  out << ">\n";
}


//namespace Impl {
//
//template <class T>
//std::uint64_t writeValuesToBuffer (std::size_t max_num_values, unsigned char* buffer,
//                                   std::vector<T> const& vec, std::size_t shift)
//{
//  std::size_t num_values = std::min(max_num_values, vec.size()-shift);
//  std::uint64_t bs = num_values*sizeof(T);
//  std::memcpy(buffer, (unsigned char*)(vec.data()+shift), std::size_t(bs));
//  return bs;
//}
//
//inline std::uint64_t compressBuffer_zlib ([[maybe_unused]] unsigned char const* buffer,
//                                          [[maybe_unused]] unsigned char* buffer_out,
//                                          [[maybe_unused]] std::uint64_t bs,
//                                          [[maybe_unused]] std::uint64_t cbs,
//                                          [[maybe_unused]] int level)
//{
//#if HAVE_VTK_ZLIB
//  uLongf uncompressed_space = uLongf(bs);
//  uLongf compressed_space = uLongf(cbs);
//
//  Bytef* out = reinterpret_cast<Bytef*>(buffer_out);
//  Bytef const* in = reinterpret_cast<Bytef const*>(buffer);
//
//  if (compress2(out, &compressed_space, in, uncompressed_space, level) != Z_OK) {
//    std::cerr << "Zlib error while compressing data.\n";
//    std::abort();
//  }
//
//  return compressed_space;
//#else
//  std::cerr << "Can not call writeCompressed without compression enabled!\n";
//  std::abort();
//  return 0;
//#endif
//}
//
//inline std::uint64_t compressBuffer_lz4 (unsigned char const* /* buffer */, unsigned char* /* buffer_out */,
//                                         std::uint64_t /* bs */, std::uint64_t /* cbs */, int /* level */)
//{
//#if HAVE_VTK_LZ4
//  std::cerr << "LZ4 Compression not yet implemented" << std::endl;
//  std::abort();
//  return 0;
//#else
//  std::cerr << "LZ4 Compression not supported. Provide the LZ4 package to CMake." << std::endl;
//  std::abort();
//  return 0;
//#endif
//}
//
//inline std::uint64_t compressBuffer_lzma (unsigned char const* /* buffer */, unsigned char* /* buffer_out */,
//                                          std::uint64_t /* bs */, std::uint64_t /* cbs */, int /* level */)
//{
//#if HAVE_VTK_LZMA
//  std::cerr << "LZMA Compression not yet implemented" << std::endl;
//  std::abort();
//  return 0;
//#else
//  std::cerr << "LZMA Compression not supported. Provide the LZMA package to CMake." << std::endl;
//  std::abort();
//  return 0;
//#endif
//}
//
//} // end namespace Impl

template <class GV, class DC>
template <class HeaderType, class FloatType>
std::uint64_t VtkWriterInterfaceMod<GV,DC>
::writeValuesAppended (std::ofstream& out, std::vector<FloatType> const& values) const
{
  assert(is_a(format_, Vtk::FormatTypes::APPENDED) && "Function should by called only in appended mode!\n");
  pos_type begin_pos = out.tellp();

  HeaderType size = values.size() * sizeof(FloatType);

  HeaderType num_full_blocks = size / block_size;
  HeaderType last_block_size = size % block_size;
  HeaderType num_blocks = num_full_blocks + (last_block_size > 0 ? 1 : 0);

  // write block-size(s)
  HeaderType zero = 0;
  if (compressor_ != Vtk::CompressorTypes::NONE) {
    out.write((char*)&num_blocks, sizeof(HeaderType));
    out.write((char*)&block_size, sizeof(HeaderType));
    out.write((char*)&last_block_size, sizeof(HeaderType));
    for (HeaderType i = 0; i < num_blocks; ++i)
      out.write((char*)&zero, sizeof(HeaderType));
  } else {
    out.write((char*)&size, sizeof(HeaderType));
  }

  HeaderType compressed_block_size = block_size + (block_size + 999)/1000 + 12;
  std::vector<unsigned char> buffer(block_size);
  std::vector<unsigned char> buffer_out;
  if (compressor_ != Vtk::CompressorTypes::NONE)
    buffer_out.resize(std::size_t(compressed_block_size));

  std::size_t num_values = block_size / sizeof(FloatType);

  std::vector<HeaderType> cbs(std::size_t(num_blocks), 0); // compressed block sizes
  for (std::size_t i = 0; i < std::size_t(num_blocks); ++i) {
    HeaderType bs = Impl::writeValuesToBuffer<FloatType>(num_values, buffer.data(), values, i*num_values);

    switch (compressor_) {
      case Vtk::CompressorTypes::NONE:
        out.write((char*)buffer.data(), bs);
        break;
      case Vtk::CompressorTypes::ZLIB:
        cbs[i] = Impl::compressBuffer_zlib(buffer.data(), buffer_out.data(), bs, compressed_block_size, compression_level);
        out.write((char*)buffer_out.data(), cbs[i]);
        break;
      case Vtk::CompressorTypes::LZ4:
        cbs[i] = Impl::compressBuffer_zlib(buffer.data(), buffer_out.data(), bs, compressed_block_size, compression_level);
        out.write((char*)buffer_out.data(), cbs[i]);
        break;
      case Vtk::CompressorTypes::LZMA:
        cbs[i] = Impl::compressBuffer_zlib(buffer.data(), buffer_out.data(), bs, compressed_block_size, compression_level);
        out.write((char*)buffer_out.data(), cbs[i]);
        break;
    }
  }

  pos_type end_pos = out.tellp();
  if (compressor_ != Vtk::CompressorTypes::NONE) {
    out.seekp(begin_pos + std::streamoff(3*sizeof(HeaderType)));
    out.write((char*)cbs.data(), std::streamsize(num_blocks*sizeof(HeaderType)));
    out.seekp(end_pos);
  }

  return std::uint64_t(end_pos - begin_pos);
}


template <class GV, class DC>
std::string VtkWriterInterfaceMod<GV,DC>
::getNames (std::vector<VtkFunction> const& data) const
{
  auto scalar = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.rangeType() == Vtk::RangeTypes::SCALAR; });
  auto vector = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.rangeType() == Vtk::RangeTypes::VECTOR; });
  auto tensor = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.rangeType() == Vtk::RangeTypes::TENSOR; });
  return (scalar != data.end() ? " Scalars=\"" + scalar->name() + "\"" : "")
      + (vector != data.end() ? " Vectors=\"" + vector->name() + "\"" : "")
      + (tensor != data.end() ? " Tensors=\"" + tensor->name() + "\"" : "");
}

} // end namespace Dune

#include <array>
#include <iosfwd>
#include <map>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/datacollectors/continuousdatacollector.hh>
#include <dune/vtk/utility/concepts.hh>

#include <dune/vtk/vtkwriterinterface.hh>

namespace Dune
{
/// File-Writer for VTK .vtu files
/**
 * Requirement:
 * - DataCollector must be a model of \ref DataCollector
 **/
template <class GridView, class DataCollector = Vtk::ContinuousDataCollector<GridView>>
class VtkUnstructuredGridWriterMod
    : public VtkWriterInterfaceMod<GridView, DataCollector>
{
  template <class> friend class VtkTimeseriesWriter;

  using Super = VtkWriterInterfaceMod<GridView, DataCollector>;
  using pos_type = typename Super::pos_type;

 public:
  /// forwarding constructor to \ref VtkWriterInterfaceMod
  using Super::Super;

 private:
  /// Write a serial VTK file in Unstructured format
  virtual void writeSerialFile (std::ofstream& out) const override;

  /// Write a parallel VTK file `pfilename.pvtu` in Unstructured format,
  /// with `size` the number of pieces and serial files given by `pfilename_p[i].vtu`
  /// for [i] in [0,...,size).
  virtual void writeParallelFile (std::ofstream& out, std::string const& pfilename, int size) const override;

  /// Write a series of timesteps in one file
  /**
   * \param filename      The name of the output file
   * \param filenameMesh  The name of a file where the mesh is stored. Must exist.
   * \param timesteps     A vector of pairs (timestep, filename) where the filename indicates
   *                      a file where the data of the timestep is stored.
   * \param blocks        A list of block sizes of the binary data stored in the files.
   *                      Order: (points, cells, pointdata[0], celldata[0], pointdata[1], celldata[1],...)
   **/
  void writeTimeseriesSerialFile (std::ofstream& out,
                                  std::string const& filenameMesh,
                                  std::vector<std::pair<double, std::string>> const& timesteps,
                                  std::vector<std::uint64_t> const& blocks) const;

  /// Write parallel VTK file for series of timesteps
  void writeTimeseriesParallelFile (std::ofstream& out,
                                    std::string const& pfilename, int size,
                                    std::vector<std::pair<double, std::string>> const& timesteps) const;

  virtual std::string fileExtension () const override
  {
    return "vtu";
  }

  virtual void writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const override;

  // Write the element connectivity to the output stream `out`. In case
  // of binary format, stores the streampos of XML attributes "offset" in the
  // vector `offsets`.
  void writeCells (std::ofstream& out,
                   std::vector<pos_type>& offsets,
                   std::optional<std::size_t> timestep = {}) const;

  void writePointIds (std::ofstream& out,
                      std::vector<pos_type>& offsets,
                      std::optional<std::size_t> timestep = {}) const;

 private:
  using Super::dataCollector_;
  using Super::format_;
  using Super::datatype_;
  using Super::headertype_;

  // attached data
  using Super::pointData_;
  using Super::cellData_;
};

// deduction guides
template <class GridView, class... Args,
    Vtk::IsGridView<GridView> = true>
VtkUnstructuredGridWriterMod(GridView, Args...)
-> VtkUnstructuredGridWriterMod<GridView, Vtk::ContinuousDataCollector<GridView>>;

template <class DataCollector, class... Args,
    Vtk::IsDataCollector<DataCollector> = true>
VtkUnstructuredGridWriterMod(DataCollector&, Args...)
-> VtkUnstructuredGridWriterMod<typename DataCollector::GridView, DataCollector>;

template <class DataCollector, class... Args,
    Vtk::IsDataCollector<DataCollector> = true>
VtkUnstructuredGridWriterMod(std::shared_ptr<DataCollector>, Args...)
-> VtkUnstructuredGridWriterMod<typename DataCollector::GridView, DataCollector>;

} // end namespace Dune


#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class GV, class DC>
void VtkUnstructuredGridWriterMod<GV,DC>
::writeSerialFile (std::ofstream& out) const
{
  std::vector<pos_type> offsets; // pos => offset
  this->writeHeader(out, "UnstructuredGrid");
  out << "<UnstructuredGrid>\n";

  out << "<Piece"
      << " NumberOfPoints=\"" << dataCollector_->numPoints() << "\""
      << " NumberOfCells=\"" << dataCollector_->numCells() << "\""
      << ">\n";

  // Write point coordinates
  out << "<Points>\n";
  this->writePoints(out, offsets);
  out << "</Points>\n";

  // Write element connectivity, types and offsets
  out << "<Cells>\n";
  writeCells(out, offsets);
  writePointIds(out, offsets);
  out << "</Cells>\n";

  // Write data associated with grid points
  out << "<PointData" << this->getNames(pointData_) << ">\n";
  for (auto const& v : pointData_)
    this->writeData(out, offsets, v, Super::POINT_DATA);
  out << "</PointData>\n";

  // Write data associated with grid cells
  out << "<CellData" << this->getNames(cellData_) << ">\n";
  for (auto const& v : cellData_)
    this->writeData(out, offsets, v, Super::CELL_DATA);
  out << "</CellData>\n";

  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";

  this->writeAppended(out, offsets);
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkUnstructuredGridWriterMod<GV,DC>
::writeParallelFile (std::ofstream& out, std::string const& pfilename, int size) const
{
  this->writeHeader(out, "PUnstructuredGrid");
  out << "<PUnstructuredGrid GhostLevel=\"0\">\n";

  // Write points
  out << "<PPoints>\n";
  out << "<PDataArray"
      << " type=\"" << to_string(datatype_) << "\""
      << " NumberOfComponents=\"3\""
      << " />\n";
  out << "</PPoints>\n";

  // Write data associated with grid points
  out << "<PPointData" << this->getNames(pointData_) << ">\n";
  for (auto const& v : pointData_) {
    out << "<PDataArray"
        << " Name=\"" << v.name() << "\""
        << " type=\"" << to_string(v.dataType()) << "\""
        << " NumberOfComponents=\"" << v.numComponents() << "\""
        << " />\n";
  }
  out << "</PPointData>\n";

  // Write data associated with grid cells
  out << "<PCellData" << this->getNames(cellData_) << ">\n";
  for (auto const& v : cellData_) {
    out << "<PDataArray"
        << " Name=\"" << v.name() << "\""
        << " type=\"" << to_string(v.dataType()) << "\""
        << " NumberOfComponents=\"" << v.numComponents() << "\""
        << " />\n";
  }
  out << "</PCellData>\n";

  // Write piece file references
  for (int p = 0; p < size; ++p) {
    std::string piece_source = pfilename + "_p" + std::to_string(p) + "." + this->fileExtension();
    out << "<Piece Source=\"" << piece_source << "\" />\n";
  }

  out << "</PUnstructuredGrid>\n";
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkUnstructuredGridWriterMod<GV,DC>
::writeTimeseriesSerialFile (std::ofstream& out,
                             std::string const& filenameMesh,
                             std::vector<std::pair<double, std::string>> const& timesteps,
                             std::vector<std::uint64_t> const& blocks) const
{
  assert(is_a(format_, Vtk::FormatTypes::APPENDED));

  std::vector<std::vector<pos_type>> offsets(timesteps.size()); // pos => offset
  this->writeHeader(out, "UnstructuredGrid");
  out << "<UnstructuredGrid"
      << " TimeValues=\"";
  {
    std::size_t i = 0;
    for (auto const& timestep : timesteps)
      out << timestep.first << (++i % 6 != 0 ? ' ' : '\n');
  }
  out << "\">\n";

  out << "<Piece"
      << " NumberOfPoints=\"" << dataCollector_->numPoints() << "\""
      << " NumberOfCells=\"" << dataCollector_->numCells() << "\""
      << ">\n";

  // Write point coordinates
  out << "<Points>\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    this->writePoints(out, offsets[i], i);
  }
  out << "</Points>\n";

  // Write element connectivity, types and offsets
  out << "<Cells>\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    writeCells(out, offsets[i], i);
    writePointIds(out, offsets[i], i);
  }
  out << "</Cells>\n";

  const std::size_t shift = offsets[0].size(); // number of blocks to write the grid

  // Write data associated with grid points
  out << "<PointData" << this->getNames(pointData_) << ">\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    for (auto const& v : pointData_)
      this->writeData(out, offsets[i], v, Super::POINT_DATA, i);
  }
  out << "</PointData>\n";

  // Write data associated with grid cells
  out << "<CellData" << this->getNames(cellData_) << ">\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    for (auto const& v : cellData_)
      this->writeData(out, offsets[i], v, Super::CELL_DATA, i);
  }
  out << "</CellData>\n";

  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";

  out << "<AppendedData encoding=\"raw\">\n_";
  pos_type appended_pos = out.tellp();

  { // write grid (points, cells)
    std::ifstream file_mesh(filenameMesh, std::ios_base::in | std::ios_base::binary);
    out << file_mesh.rdbuf();
    assert( std::uint64_t(out.tellp()) == std::accumulate(blocks.begin(), std::next(blocks.begin(),shift), std::uint64_t(appended_pos)) );
  }

  // write point-data and cell-data
  for (auto const& timestep : timesteps) {
    std::ifstream file(timestep.second, std::ios_base::in | std::ios_base::binary);
    out << file.rdbuf();
  }
  out << "</AppendedData>\n";

  out << "</VTKFile>";

  // write correct offsets in file.
  pos_type offset = 0;
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    offset = 0;
    auto const& off = offsets[i];

    // write mesh data offsets
    for (std::size_t j = 0; j < shift; ++j) {
      out.seekp(off[j]);
      out << '"' << offset << '"';
      offset += pos_type(blocks[j]);
    }
  }

  std::size_t j = shift;
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    auto const& off = offsets[i];

    for (std::size_t k = shift; k < off.size(); ++k) {
      out.seekp(off[k]);
      out << '"' << offset << '"';
      offset += pos_type(blocks[j++]);
    }
  }
}


template <class GV, class DC>
void VtkUnstructuredGridWriterMod<GV,DC>
::writeTimeseriesParallelFile (std::ofstream& out,
                               std::string const& pfilename,
                               int size,
                               std::vector<std::pair<double, std::string>> const& timesteps) const
{
  this->writeHeader(out, "PUnstructuredGrid");
  out << "<PUnstructuredGrid GhostLevel=\"0\""
      << " TimeValues=\"";
  {
    std::size_t i = 0;
    for (auto const& timestep : timesteps)
      out << timestep.first << (++i % 6 != 0 ? ' ' : '\n');
  }
  out << "\">\n";

  // Write points
  out << "<PPoints>\n";
  out << "<PDataArray"
      << " type=\"" << to_string(datatype_) << "\""
      << " NumberOfComponents=\"3\""
      << " />\n";
  out << "</PPoints>\n";

  // Write data associated with grid points
  out << "<PPointData" << this->getNames(pointData_) << ">\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    for (auto const& v : pointData_) {
      out << "<PDataArray"
          << " Name=\"" << v.name() << "\""
          << " type=\"" << to_string(v.dataType()) << "\""
          << " NumberOfComponents=\"" << v.numComponents() << "\""
          << " TimeStep=\"" << i << "\""
          << " />\n";
    }
  }
  out << "</PPointData>\n";

  // Write data associated with grid cells
  out << "<PCellData" << this->getNames(cellData_) << ">\n";
  for (std::size_t i = 0; i < timesteps.size(); ++i) {
    for (auto const& v : cellData_) {
      out << "<PDataArray"
          << " Name=\"" << v.name() << "\""
          << " type=\"" << to_string(v.dataType()) << "\""
          << " NumberOfComponents=\"" << v.numComponents() << "\""
          << " TimeStep=\"" << i << "\""
          << " />\n";
    }
  }
  out << "</PCellData>\n";

  // Write piece file references
  for (int p = 0; p < size; ++p) {
    std::string piece_source = pfilename + "_p" + std::to_string(p) + "." + this->fileExtension();
    out << "<Piece Source=\"" << piece_source << "\" />\n";
  }

  out << "</PUnstructuredGrid>\n";
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkUnstructuredGridWriterMod<GV,DC>
::writeCells (std::ofstream& out, std::vector<pos_type>& offsets,
              std::optional<std::size_t> timestep) const
{
  if (format_ == Vtk::FormatTypes::ASCII) {
    auto cells = dataCollector_->cells();
    out << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\"";
    if (timestep)
      out << " TimeStep=\"" << *timestep << "\"";
    out << ">\n";
    this->writeValuesAscii(out, cells.connectivity);
    out << "</DataArray>\n";

    out << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\"";
    if (timestep)
      out << " TimeStep=\"" << *timestep << "\"";
    out << ">\n";
    this->writeValuesAscii(out, cells.offsets);
    out << "</DataArray>\n";

    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"";
    if (timestep)
      out << " TimeStep=\"" << *timestep << "\"";
    out << ">\n";
    this->writeValuesAscii(out, cells.types);
    out << "</DataArray>\n";
  }
  else { // Vtk::FormatTypes::APPENDED format
    out << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\"";
    if (timestep)
      out << " TimeStep=\"" << *timestep << "\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\"";
    if (timestep)
      out << " TimeStep=\"" << *timestep << "\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\"";
    if (timestep)
      out << " TimeStep=\"" << *timestep << "\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
void VtkUnstructuredGridWriterMod<GV,DC>
::writePointIds (std::ofstream& out,
                 std::vector<pos_type>& offsets,
                 std::optional<std::size_t> timestep) const
{
  auto ids = dataCollector_->pointIds();
  if (ids.empty())
    return;

  if (format_ == Vtk::FormatTypes::ASCII) {
    out << "<DataArray type=\"UInt64\" Name=\"global_point_ids\" format=\"ascii\"";
    if (timestep)
      out << " TimeStep=\"" << *timestep << "\"";
    out << ">\n";
    this->writeValuesAscii(out, ids);
    out << "</DataArray>\n";
  }
  else { // Vtk::FormatTypes::APPENDED format
    out << "<DataArray type=\"UInt64\" Name=\"global_point_ids\" format=\"appended\"";
    if (timestep)
      out << " TimeStep=\"" << *timestep << "\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
void VtkUnstructuredGridWriterMod<GV,DC>
::writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const
{
  assert(is_a(format_, Vtk::FormatTypes::APPENDED) && "Function should by called only in appended mode!\n");

  Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(datatype_, headertype_,
                                                              [&](auto f, auto h) {
                                                                using F = typename decltype(f)::type;
                                                                using H = typename decltype(h)::type;

                                                                // write points
                                                                blocks.push_back(this->template writeValuesAppended<H>(out, dataCollector_->template points<F>()));

                                                                // write connectivity, offsets, and types
                                                                auto cells = dataCollector_->cells();
                                                                blocks.push_back(this->template writeValuesAppended<H>(out, cells.connectivity));
                                                                blocks.push_back(this->template writeValuesAppended<H>(out, cells.offsets));
                                                                blocks.push_back(this->template writeValuesAppended<H>(out, cells.types));

                                                                // optionally, write global point IDs
                                                                auto ids = dataCollector_->pointIds();
                                                                if (!ids.empty())
                                                                  blocks.push_back(this->template writeValuesAppended<H>(out, ids));
                                                              });
}

} // end namespace Dune

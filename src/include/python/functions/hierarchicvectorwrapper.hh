// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_FUNCTIONS_HIERARCHICVECTORWRAPPER_HH
#define DUNE_PYTHON_FUNCTIONS_HIERARCHICVECTORWRAPPER_HH

#include <type_traits>
#include <utility>

#include <dune/typetree/utility.hh>

#include <dune/functions/common/indexaccess.hh>

namespace Dune
{

  namespace Python
  {

    // HierarchicalVectorWrapper
    // -------------------------

    template< class V, class C, class Holder = V >
    class HierarchicVectorWrapper
    {
      typedef HierarchicVectorWrapper< V, C, Holder > This;

    public:
      typedef V Vector;

      typedef typename Vector::size_type size_type;

      template< class MultiIndex >
      using Entry = C;

      template< class... Args, std::enable_if_t< std::is_constructible< Holder, Args &&... >::value, int > = 0 >
      HierarchicVectorWrapper ( Args &&... args )
        : holder_( std::forward< Args >( args )... )
      {}

      template< class MultiIndex >
      const Entry< MultiIndex > &operator[] ( const MultiIndex &index ) const
      {
        return Dune::Functions::hybridMultiIndexAccess< const Entry< MultiIndex > & >( vector(), index );
      }

      template< class MultiIndex >
      Entry< MultiIndex > &operator[] ( const MultiIndex &index )
      {
        return Dune::Functions::hybridMultiIndexAccess< Entry< MultiIndex > & >( vector(), index );
      }

      template< class MultiIndex >
      const Entry< MultiIndex > &operator() ( const MultiIndex &index ) const
      {
        return (*this)[ index ];
      }

      template< class MultiIndex >
      Entry< MultiIndex > &operator() ( const MultiIndex &index )
      {
        return (*this)[ index ];
      }

      const Vector &vector () const { return get( holder_ ); }
      Vector &vector () { return get( holder_ ); }

    private:
      template< class H >
      static decltype( static_cast< const Vector & >( std::declval< const H & >() ) ) get ( const H &holder )
      {
        return static_cast< const Vector & >( holder );
      }

      template< class H >
      static decltype( static_cast< Vector & >( std::declval< H & >() ) ) get ( H &holder )
      {
        return static_cast< Vector & >( holder );
      }

      template< class H >
      static decltype( static_cast< const Vector & >( *std::declval< const H & >() ) ) get ( const H &holder )
      {
        return static_cast< const Vector & >( *holder );
      }

      template< class H >
      static decltype( static_cast< Vector & >( *std::declval< H & >() ) ) get ( H &holder )
      {
        return static_cast< Vector & >( *holder );
      }

      Holder holder_;
    };

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_FUNCTIONS_HIERARCHICVECTORWRAPPER_HH

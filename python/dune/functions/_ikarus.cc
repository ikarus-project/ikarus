#include <config.h>

#include <cmath>
#include <sstream>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/python/pybind11/extensions.h>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( _ikarus, module )
{
  try
  {
    {
      typedef typename Dune::Fem::detail::SingletonStorage::StorageType Singleton;
      pybind11::class_< Singleton, std::shared_ptr<Singleton> > cls( module, "_Singleton" );
      Dune::Fem::detail::SingletonStorage::getStorage();
      module.attr( "_singleton" ) = pybind11::cast( Dune::Fem::detail::SingletonStorage::storage_ );
    }
    int argc = 0;
    char **argv = nullptr;
    Dune::Fem::MPIManager::initialize( argc, argv );

    if( !pybind11::already_registered< Dune::Fem::MPIManager::Communication >() )
      DUNE_THROW( Dune::Exception, "Communication not registered, yet" );

    module.attr( "comm" ) = pybind11::cast( Dune::Fem::MPIManager::comm() );
  }
  catch ( const std::exception &e )
  {
    std::cout << e.what() << std::endl;
  }

  {
    using pybind11::operator""_a;
    using pybind11::str;

    auto mpiManagerCls = pybind11::class_<Dune::Fem::MPIManager>(module, "threading");
    pybind11::class_< Dune::Fem::ParameterContainer > param( module, "Parameter" );

    // function to verbose rank and verbosity level
    param.def( "_setVerbosity", [] ( Dune::Fem::ParameterContainer &self, const int level, const int rank )
        {
          // set verbosity level (see dune/fem/io/parameter/container.hh for doc)
          self.append( "fem.verbositylevel", level, /* force */ true );
          // set verbose rank
          self.append( "fem.verboserank", rank, /* force */ true );
        }, "level"_a, "rank"_a );

    param.def( "write", [] ( Dune::Fem::ParameterContainer &self, const std::string &fileName ) {
          std::ofstream file( fileName );
          if( file )
          {
            self.write( file );
            file.close();
          }
        }, "fileName"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &fileName ) {
          self.append( fileName );
        }, "fileName"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::map< std::string, std::string > &entries ) {
          for( auto entry : entries )
            self.append( entry.first, entry.second );
        }, "entries"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const pybind11::dict &entries ) {
          for ( auto entry : entries )
          {
            std::string s = str( entry.second );
            if (s == "False")
              self.append( str( entry.first ), "false" );
            else if (s == "True")
              self.append( str( entry.first ), "true" );
            else
              self.append( str( entry.first ), s );
          }
        }, "entries"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, pybind11::handle value ) {
          self.append( key, str( value ) );
        }, "key"_a, "value"_a );

#if 0
    // do we really need this one?
    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, int value ) {
        self.append( key, std::to_string( value ) );
      }, "key"_a, "value"_a );

    // do we really need this one?
    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, double value ) {
        self.append( key, std::to_string( value ) );
      }, "key"_a, "value"_a );
#endif

    param.def( "exists", [] ( const Dune::Fem::ParameterContainer &self, const std::string &key ) {
          return self.exists( key );
        }, "key"_a );

    param.def( "__getitem__", [] ( const Dune::Fem::ParameterContainer &self, const std::string &key ) {
          if (!self.exists( key ))
            throw pybind11::key_error("key not found in parameter file");
          return self.getValue< std::string >( key );
        } , "key"_a );

    param.def( "__setitem__", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, pybind11::handle value ) {
          self.append( key, str( value ) );
        }, "key"_a, "value"_a );

    param.def( "get", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, bool defaultValue ) -> bool {
          if (self.exists( key ))
            return self.getValue< bool >( key );
          else
            self.append( key, defaultValue?"true":"false" );
          return defaultValue;
        }, "key"_a, "defaultValue"_a );
    param.def( "get", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, int defaultValue ) {
          if (self.exists( key ))
            return self.getValue< int >( key );
          else
            self.append( key, std::to_string(defaultValue) );
          return defaultValue;
        }, "key"_a, "defaultValue"_a );
    param.def( "get", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, double defaultValue ) {
          if (self.exists( key ))
            return self.getValue< double >( key );
          else
            self.append( key, std::to_string(defaultValue) );
          return defaultValue;
        }, "key"_a, "defaultValue"_a );

    param.def( "get", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, pybind11::handle defaultValue ) {
          if (self.exists( key ))
            return str( self.getValue< std::string >( key ) );
          else if (!defaultValue.is_none())
            self.append( key, str( defaultValue ) );
          else
            throw pybind11::key_error("key not found in parameter file");
          return str( defaultValue );
        }, "key"_a, pybind11::arg("defaultValue")=pybind11::none() );

    param.def( "__str__", [] ( const Dune::Fem::ParameterContainer &self ) {
      std::stringstream s;
      self.write( s );
      return s.str();
    } );

    module.attr( "parameter" ) = pybind11::cast( Dune::Fem::Parameter::container(),
                                              pybind11::return_value_policy::reference );

    // add finalize method for MPI and PETSc
    module.def( "__finalizeFemModule__", [] () {
      pybind11::module gc = pybind11::module::import("gc");
      gc.attr("collect")();
      Dune::Fem::MPIManager::comm().barrier();
      Dune::Fem::MPIManager::finalize();
    } );

    mpiManagerCls.def_property_readonly_static("max",
                                               [](pybind11::object){ return Dune::Fem::MPIManager::maxThreads(); });
    mpiManagerCls.def_property_static("use",
        [](pybind11::object){ return Dune::Fem::MPIManager::numThreads(); },
        [](pybind11::object,uint threads){ Dune::Fem::MPIManager::setNumThreads(threads); });
    mpiManagerCls.def_static("useMax",
                             [](){ Dune::Fem::MPIManager::setNumThreads(Dune::Fem::MPIManager::maxThreads()); });

  }
}

/* begin Ikarus
put the definitions for config.h specific to
your project here. Everything above will be
overwritten
*/

/* begin private */
/* Name of package */
#define PACKAGE "@DUNE_MOD_NAME@"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "@DUNE_MAINTAINER@"

/* Define to the full name of this package. */
#define PACKAGE_NAME "@DUNE_MOD_NAME@"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "@DUNE_MOD_NAME@ @DUNE_MOD_VERSION@"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "@DUNE_MOD_NAME@"

/* Define to the home page for this package. */
#define PACKAGE_URL "@DUNE_MOD_URL@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@DUNE_MOD_VERSION@"

/* end private */

/* Define to the version of Ikarus */
#define IKARUS_VERSION "@IKARUS_VERSION@"

/* Define the sparse matrix addon for eigen */
#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"

/* Init eigen matrices with nan */
#define EIGEN_INITIALIZE_MATRICES_BY_NAN

/* Define to the major version of Ikarus */
#define IKARUS_VERSION_MAJOR @IKARUS_VERSION_MAJOR@

/* Define to the minor version of Ikarus */
#define IKARUS_VERSION_MINOR @IKARUS_VERSION_MINOR@

/* Define to the revision of Ikarus */
#define IKARUS_VERSION_REVISION @IKARUS_VERSION_REVISION@

/* end Ikarus
Everything below here will be overwritten
*/

//------------------------------------------------------------------------------

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
// from LINALG
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"
#include "fn_r2.hpp"
#include "fn_norm.hpp"
#include "HDF5_Tools.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"


//------------------------------------------------------------------------------
// from MTK
#include "cl_MTK_Mesh.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>


//------------------------------------------------------------------------------
// HMR
#define private public
#define protected public
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Mesh.hpp"


#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#undef private
#undef protected

//------------------------------------------------------------------------------

#include "cl_MTK_Mapper.hpp"

// select namespaces
using namespace moris;
using namespace hmr;

real
CircleFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

real
SimionescuFunction( const Matrix< DDRMat > & aPoint )
{
    return 0.1 * aPoint( 0 ) * aPoint ( 1 );
}

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
//------------------------------------------------------------------------------

/*!
 * \section Facets Tutorial
 * This example creates a simple 2x2 mesh and tests MTK functionality.
 * This example was written for Keenan and will be turned into a test soon.
 */
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

//------------------------------------------------------------------------------

    // get root from environment
    std::string tMorisRoot = std::getenv("MORISROOT");

    // determine path for object file
    std::string tHdf5FilePath = tMorisRoot + "/projects/HMR/test/data/hmr_sideset_test_2d.hdf5" ;

//------------------------------------------------------------------------------


      ParameterList tParameters = create_hmr_parameter_list();

      tParameters.set( "number_of_elements_per_dimension", "4, 6" );

      tParameters.set( "domain_offset", "0, 0" );
      tParameters.set( "domain_dimensions", "4, 6" );
      tParameters.set( "domain_sidesets", "1, 2, 3, 4" );
      tParameters.set( "bspline_orders", "1" );
      tParameters.set( "lagrange_orders", "1" );
      tParameters.set( "verbose", 0 );

//------------------------------------------------------------------------------

      HMR tHMR( tParameters );

//------------------------------------------------------------------------------
//    create refinement pattern
//------------------------------------------------------------------------------
      std::shared_ptr< Database > tDatabase = tHMR.get_database();

      tDatabase->set_activation_pattern( tHMR.get_parameters()->get_output_pattern() );

      for( uint tLevel = 0; tLevel < 4; ++tLevel )
      {
          // flag first element
          tDatabase->flag_element( 0 );

          // flag last element
          tDatabase->flag_element(
                  tDatabase->get_number_of_elements_on_proc()-1 );

          // manually refine, do not reset pattern
          tDatabase->get_background_mesh()->perform_refinement();

      }

      // update database etc
      tDatabase->perform_refinement( false );
      //tHMR.perform_refinement();

      // finish mesh
      tHMR.finalize();
//------------------------------------------------------------------------------
//    create Mesh and variables
//------------------------------------------------------------------------------
      // create MTK mesh
      std::shared_ptr< Mesh > tMesh = tHMR.create_mesh();

      Matrix<IndexMat> tElements;
      Matrix<IndexMat> tElementsSolution;
      Matrix<IndexMat> tOrds;

//------------------------------------------------------------------------------
//    write solution
//------------------------------------------------------------------------------

     // create file
     hid_t tFileID = create_hdf5_file( tHdf5FilePath );

      // error handler
      herr_t tStatus;

      for( uint s=1; s<=4; ++s )
      {
          // create label
          std::string tSetLabel = "SideSet_" + std::to_string( s );

          // ask mesh for sideset
          tMesh->get_sideset_elems_loc_inds_and_ords(
                  tSetLabel, tElements, tOrds );

          // save data
          save_matrix_to_hdf5_file( tFileID, tSetLabel, tElements, tStatus );
      }


      // close file
      close_hdf5_file( tFileID );

      // save exodus file for visual inspection
      tHMR.save_to_exodus( "Mesh.exo" );

//------------------------------------------------------------------------------
//    open solution
//------------------------------------------------------------------------------

/*      // create file
      hid_t tFileID = open_hdf5_file( tHdf5FilePath );

      // error handler
      herr_t tStatus = 0;

      for( uint s=1; s<=6; ++s )
      {
          // create label
          std::string tSetLabel = "SideSet_" + std::to_string( s );

          // ask mesh for sideset
          tMesh->get_sideset_elems_loc_inds_and_ords(
                  tSetLabel, tElements, tOrds );

          // read
          load_matrix_from_hdf5_file( tFileID, tSetLabel, tElementsSolution, tStatus );

          if( tElementsSolution.length() > 0 )
          {
              // compare result
              //REQUIRE( all_true( tElementsAtSurface == tElementsAtSurfaceExpect ) );
              std::cout << "REQUIRE " <<  all_true( tElementsSolution == tElements ) << std::endl;
          }
      }

      // close file
      close_hdf5_file( tFileID ); */

//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

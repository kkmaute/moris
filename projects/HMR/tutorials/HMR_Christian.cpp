//------------------------------------------------------------------------------

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_MTK_Mesh.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------
// LINALG

#include "fn_r2.hpp"
#include "fn_norm.hpp"

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


    ParameterList tParameters = create_hmr_parameter_list();

      tParameters.set( "number_of_elements_per_dimension", "5, 3, 2" );

      tParameters.set( "domain_offset", "-2, -2, -2" );
      tParameters.set( "domain_dimensions", "5, 3, 2" );

      tParameters.set( "bspline_orders", "2" );
      tParameters.set( "lagrange_orders", "2" );
      tParameters.set( "verbose", 1 );

//------------------------------------------------------------------------------

      HMR tHMR( tParameters );

//------------------------------------------------------------------------------

      tHMR.finalize();

      Cell< Background_Element_Base* > tElements;
      tHMR.mDatabase->get_background_mesh()->collect_coarsest_elements_on_side(
              2,
              tElements
              );

      for( Background_Element_Base* tElement : tElements )
      {
          std::cout << tElement->get_domain_id() << std::endl;
      }

//------------------------------------------------------------------------------

/*      // create mesh
      auto tMesh = tHMR.create_mesh();

      uint tOrder = 2;

      std::shared_ptr< Field > tField = tMesh->create_field( "Circle", tOrder );
      tField->evaluate_scalar_function( CircleFunction );

//------------------------------------------------------------------------------

      // create mapper with one mesh
      mapper::Mapper tMapper( tMesh.get() ); // < -- also add two meshes if desired

      // map node to B-Splines
      tMapper.perform_mapping(
              "Circle",
              EntityRank::NODE,
              "Circle",
              EntityRank::BSPLINE_2 );

      // map B-Splines to Nodes
      tMapper.perform_mapping(
                    "Circle",
                    EntityRank::BSPLINE_2,
                    "Circle",
                    EntityRank::NODE );

//------------------------------------------------------------------------------

      tField->save_field_to_hdf5("Circle.hdf5");  */

      tHMR.save_background_mesh_to_vtk("BG.vtk");
      tHMR.save_to_exodus( 1, "Mesh.exo" );

//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

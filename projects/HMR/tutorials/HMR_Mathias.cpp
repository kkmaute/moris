//------------------------------------------------------------------------------
#include <memory>
#include <string>

// dynamik linker function
#include "dlfcn.h"

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
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_Mesh_Factory.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>


//------------------------------------------------------------------------------
// HMR

#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Mesh.hpp"


#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"




//------------------------------------------------------------------------------

// select namespaces
using namespace moris;
using namespace hmr;


//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
//------------------------------------------------------------------------------

real
CircleFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}
//------------------------------------------------------------------------------
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

//------------------------------------------------------------------------------
// this example creates a mesh and performs a manual refinement
//------------------------------------------------------------------------------

    // order for this example
    uint tOrder = 1;

    // create parameter object
    Parameters tParameters;
    tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );
    tParameters.set_verbose( false );
    tParameters.set_multigrid( true );
    tParameters.set_mesh_orders_simple( tOrder );

    // create HMR object
    HMR tHMR( tParameters );

    // flag first element for refinement
    tHMR.flag_element( 0 );
    tHMR.perform_refinement( gRefinementModeBSpline );

    tHMR.flag_element( 0 );
    tHMR.perform_refinement( gRefinementModeBSpline );

    tHMR.finalize();

    // grab the output mesh
    std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tOrder );

    // create a field on this mesh
    std::shared_ptr< Field > tField = tMesh->create_field( "Circle", tOrder );
    tField->evaluate_scalar_function( CircleFunction );


    std::cout << std::endl;
    // List of B-Spline Meshes
    for( uint k=0; k<tHMR.get_database()->get_number_of_bspline_meshes(); ++k )
    {
        // get pointer to mesh
        BSpline_Mesh_Base * tMesh = tHMR.get_database()->get_bspline_mesh_by_index( k );

        std::cout << "BSpline Mesh " << k <<
                ": active pattern " << tMesh->get_activation_pattern() <<
                " order " << tMesh->get_order() <<
                " active basis " << tMesh->get_number_of_active_basis_on_proc()
                << std::endl;

    }
    std::cout << std::endl;

    // List of Lagrange Meshes
    for( uint k=0; k<tHMR.get_database()->get_number_of_lagrange_meshes(); ++k )
    {
        // get pointer to mesh
        Lagrange_Mesh_Base * tMesh = tHMR.get_database()->get_lagrange_mesh_by_index( k );

        std::cout << "Lagrange Mesh " << k <<
                ": active pattern " << tMesh->get_activation_pattern() <<
                " order " << tMesh->get_order() <<
                " active basis " << tMesh->get_number_of_nodes_on_proc() << std::endl;

    }
    std::cout << std::endl;

    map< moris_id, moris_index > tMap;
    tMesh->get_adof_map( tOrder, tMap );
    tMap.print("Adof Map");

    tHMR.save_bsplines_to_vtk("BSplines.vtk");
    tHMR.save_to_exodus("Mesh.exo");
//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

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
#define protected public
#define private   public
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#undef protected
#undef private
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
    return ( norm( aPoint ) - 1.2 );
}

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );
//------------------------------------------------------------------------------

    ParameterList tParameters = create_hmr_parameter_list();

    tParameters.set( "number_of_elements_per_dimension", "2,2,2" );
    tParameters.set( "domain_dimensions", "2, 2,2" );
    tParameters.set( "domain_offset", "0, 0, 0" );
    tParameters.set( "verbose", 1 );
    tParameters.set( "truncate_bsplines", 1 );

    tParameters.set( "bspline_orders", "1" );
    tParameters.set( "lagrange_orders", "1" );



    //HMR tHMR( "pattern_1.hdf5", "pattern_2.hdf5" );

    HMR tHMRA( tParameters );

    for( uint tLevel = 0; tLevel < 2; ++tLevel )
    {
        tHMRA.get_database()->flag_element( 0 );

        // manually refine, do not reset pattern
        tHMRA.perform_refinement();
    }

    tHMRA.save_to_hdf5( "pattern_1.hdf5" ) ;

    HMR tHMRB( tParameters );

    for( uint tLevel = 0; tLevel < 2; ++tLevel )
    {
        tHMRB.get_database()->flag_element( tHMRB.get_database()->get_number_of_elements_on_proc()-1  );

        // manually refine, do not reset pattern
        tHMRB.perform_refinement();
    }

    tHMRB.save_to_hdf5( "pattern_2.hdf5" );

    HMR tHMR( "pattern_1.hdf5", "pattern_2.hdf5" );

    tHMR.save_background_mesh_to_vtk( "Background.vtk" );
    tHMR.save_faces_to_vtk( "Faces.vtk" );
//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;
}

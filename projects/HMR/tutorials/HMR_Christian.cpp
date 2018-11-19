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

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

//------------------------------------------------------------------------------

    ParameterList tParameters = create_hmr_parameter_list();

    tParameters.set( "number_of_elements_per_dimension", "4, 4" );
    tParameters.set( "domain_dimensions", "2, 2" );
    tParameters.set( "domain_offset", "-1, -1" );

    tParameters.set( "verbose", 1 );
    tParameters.set( "truncate_bsplines", 1 );

    tParameters.set( "bspline_orders", "1" );
    tParameters.set( "lagrange_orders", "1" );

    // create HMR object
    HMR tHMRA( tParameters );

    tHMRA.flag_element( 0 );
    tHMRA.perform_refinement( true );

    tHMRA.flag_element( 0 );
    tHMRA.perform_refinement( false );

    tHMRA.finalize();
    tHMRA.save_mesh_to_vtk("Lagrange.vtk");
    tHMRA.save_bsplines_to_vtk("BSpline.vtk");

    tHMRA.get_database()->set_activation_pattern( 0 );
    tHMRA.save_to_hdf5( "MyDatabase.hdf5");

    std::cout << "----------------" << std::endl;

    HMR tHMRB( "MyDatabase.hdf5" );
    tHMRB.finalize();

    tHMRB.save_mesh_to_vtk("Lagrange2.vtk");
    tHMRB.save_bsplines_to_vtk("BSpline2.vtk");

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

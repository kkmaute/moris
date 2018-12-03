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
#include "fn_sum.hpp"

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

#include "cl_HMR_Paramfile.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Mesh.hpp"


#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_BSpline.hpp"

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

<<<<<<< HEAD
//------------------------------------------------------------------------------

    std::string tPath = "/home/messe/Examples/HMR_Opt_Test_2d/parameters.xml";

    Paramfile tParams( tPath, State::INITIALIZE_MESH );

=======

		//uint tMaxLevel = 1;
		uint tOrder = 2;
		uint tDimension = 3;

//------------------------------------------------------------------------------
		// create settings object
		auto tParameters = std::make_shared< Parameters >();

		// set number of elements
		moris::Matrix< moris::DDLUMat > tNumberOfElements;

		tNumberOfElements.set_size( tDimension, 1, 4 );

		tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

		// do not print debug information during test
		tParameters->set_verbose( true );

		// deactivate truncation
		tParameters->set_bspline_truncation( false );

		tParameters->set_multigrid( true );

		// set buffer size to zero
		tParameters->set_buffer_size( 0 );

		// create factory
		moris::hmr::Factory tFactory;

		// set buffer size to zero
		tParameters->set_buffer_size( tOrder );
		tParameters->set_additional_lagrange_refinement( 1 );

		// set aura
		//tParameters->set_max_polynomial( tOrder );

		// set simple mesh order
		tParameters->set_mesh_orders_simple( tOrder );

		HMR tHMR( tParameters.get() );



		tHMR.flag_element( 0 );
		tHMR.perform_initial_refinement();


		tHMR.finalize();
		tHMR.save_bsplines_to_vtk("Basis.vtk");

		std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tOrder );

		std::shared_ptr< Field > tOutputField = tMesh->create_field( "MyField", tOrder ) ;

		// interpolate field onto union mesh
		/*tHMR.get_database()->interpolate_field(
		        tHMR.get_parameters()->get_lagrange_input_pattern(),
		        tInputField,
		        tHMR.get_parameters()->get_lagrange_output_pattern(),
		        tOutputField ); */


		Matrix< IndexMat > tElements = tMesh->get_entity_connected_to_entity_loc_inds( 0,
		        EntityRank::ELEMENT,
		        EntityRank::BSPLINE_2 );

		print( tElements, "tElements" );

//------------------------------------------------------------------------------
		//delete tParameters;
>>>>>>> dcb9c5dee35f43d79480f100eb0ce0935a7973e6
//------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;

}

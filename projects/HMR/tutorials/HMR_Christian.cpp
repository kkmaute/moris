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

//real
//nchoosek( const uint & aN, const uint aK )
//{
//    real aResult = 1.0;
//
//    for ( uint i=1; i<=aK; ++i )
//    {
//        aResult *= ( ( real ) aN+1-i ) / ( real( i ) );
//    }
//
//    return aResult;
//}

//------------------------------------------------------------------------------

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );


		//uint tMaxLevel = 1;
		uint tOrder = 2;
		uint tDimension = 3;

//------------------------------------------------------------------------------
		// create settings object
		auto tParameters = std::make_shared< Parameters >();

		// set number of elements
		moris::Matrix< moris::DDLUMat > tNumberOfElements;

		tNumberOfElements.set_size( tDimension, 1, 1 );

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

		// set aura
		//tParameters->set_max_polynomial( tOrder );

		// set simple mesh order
		tParameters->set_mesh_orders_simple( tOrder );

		HMR tHMR( tParameters.get() );

		std::shared_ptr< Field > tInputField = tHMR.create_field( "MyField" );

		Matrix< DDRMat > & tNodeValues = tInputField->get_node_values();

		tNodeValues.set_size( 27, 1 );

		tNodeValues(0) = 0.9296875;
		tNodeValues(1) = 0.887709466748001;
		tNodeValues(2) = 0.890625;
		tNodeValues(3) = 0.9375;
		tNodeValues(4) = 0.71875;
		tNodeValues(5) = 0.554045894498667;
		tNodeValues(6) = 0.5625;
		tNodeValues(7) = 0.75;
		tNodeValues(8) = 0.905521116687;
		tNodeValues(9) = 0.889896116687;
		tNodeValues(10) = 0.91015625;
		tNodeValues(11) = 0.935546875;
		tNodeValues(12) = 0.841796875;
		tNodeValues(13) = 0.748629511185668;
		tNodeValues(14) = 0.75390625;
		tNodeValues(15) = 0.859375;
		tNodeValues(16) = 0.622886473624667;
		tNodeValues(17) = 0.560386473624667;
		tNodeValues(18) = 0.640625;
		tNodeValues(19) = 0.7421875;
		tNodeValues(20) = 0.795324500699104;
		tNodeValues(21) = 0.90899746667175;
		tNodeValues(22) = 0.636190368406167;
		tNodeValues(23) = 0.85498046875;
		tNodeValues(24) = 0.752587065296417;
		tNodeValues(25) = 0.787743315296417;
		tNodeValues(26) = 0.7978515625;

		tHMR.flag_element( 0 );
		tHMR.perform_refinement( RefinementMode::SIMPLE );


		tHMR.finalize();

		std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tOrder );

		std::shared_ptr< Field > tOutputField = tMesh->create_field( "MyField", tOrder ) ;

		// interpolate field onto union mesh
		tHMR.get_database()->interpolate_field(
		        tHMR.get_parameters()->get_lagrange_input_pattern(),
		        tInputField,
		        tHMR.get_parameters()->get_lagrange_output_pattern(),
		        tOutputField );

		tHMR.save_to_exodus("Mesh.exo");
		tHMR.save_to_exodus( 0, "LastStep.exo");
//------------------------------------------------------------------------------
		//delete tParameters;
//------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;

}

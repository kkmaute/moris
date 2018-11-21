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


		uint tMaxLevel = 1;
		uint tOrder = 2;
		uint tDimension = 2;

//------------------------------------------------------------------------------
		// create settings object
		moris::hmr::Parameters * tParameters = new moris::hmr::Parameters;

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

		// create background mesh object
		moris::hmr::Background_Mesh_Base* tBackgroundMesh
		= tFactory.create_background_mesh( tParameters );

		// refine a few elements in the mesh
		for( moris::uint l=0; l<tMaxLevel; ++l  )
		{
			auto tNumberOfElements
			=  tBackgroundMesh->get_number_of_active_elements_on_proc();

			// refine every other element
			for( moris::luint k=0; k<tNumberOfElements; ++k )
			{
				// get element
				moris::hmr::Background_Element_Base* tElement
					= tBackgroundMesh->get_element( k );

				// flag element for refinement
				tElement->put_on_refinement_queue();
			}

			// refine mesh
			tBackgroundMesh->perform_refinement();
		}

		// create B-Spline mesh
		moris::hmr::BSpline_Mesh_Base* tBSplineMesh
			= tFactory.create_bspline_mesh( tParameters, tBackgroundMesh, 0, tOrder );

		tBSplineMesh->update_mesh();
		Matrix< IdMat > tCommTable;
		tBSplineMesh->calculate_basis_indices( tCommTable );

		tBackgroundMesh->save_to_vtk( "BackgroundMesh.vtk" );
		tBSplineMesh->save_to_vtk("BSplines.vtk");

		Basis * tBasis = tBSplineMesh->get_basis_by_index( 4 );

		std::cout << "Basis " << tBasis->get_domain_id() << std::endl;



		uint tNumberOfChildrenPerDirection = tOrder + 2;
		uint tNumberOfChildren = std::pow( tNumberOfChildrenPerDirection, tDimension );

		Matrix< DDRMat > tWeights( tNumberOfChildren, 1 );

		uint tChild = 0;

		for ( uint j=0; j< tNumberOfChildrenPerDirection ; ++j )
		{
			for( uint i=0; i<tNumberOfChildrenPerDirection; ++i )
			{
				tWeights( tChild++ ) = nchoosek( tOrder+1, i ) * nchoosek( tOrder+1, j );
			}
		}

		tWeights = tWeights / std::pow( 2, tOrder );

		for( uint c=0; c<tNumberOfChildren; ++c )
		{
			std::cout << "   Child " << c
					  << " id:" << tBasis->get_child( c )->get_domain_id()
					  << " weight:" << tWeights( c ) << std::endl;
		}
//------------------------------------------------------------------------------
		delete tBSplineMesh;
		delete tBackgroundMesh;
		delete tParameters;

//------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;

}

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


       //moris::uint tMyCoeff = 1;

       std::cout<<"---"<<std::endl;



       ParameterList tParameters = create_hmr_parameter_list();

       tParameters.set( "number_of_elements_per_dimension", "2, 2" );
       tParameters.set( "domain_dimensions", "3, 3" );
       tParameters.set( "domain_offset", "-1.5, -1.5" );
       tParameters.set( "verbose", 1 );
       tParameters.set( "truncate_bsplines", 1 );
       tParameters.set( "bspline_orders", "2" );
       tParameters.set( "lagrange_orders", "2" );
       tParameters.set( "additional_lagrange_refinement", 2 );

       HMR tHMR( tParameters );

       // std::shared_ptr< Database >
       //auto tDatabase = tHMR.get_database();

       // manually select output pattern
       //tDatabase->get_background_mesh()->set_activation_pattern( tHMR.get_parameters()->get_lagrange_output_pattern() );

       tHMR.perform_initial_refinement();

       //tDatabase->get_background_mesh()->get_element(0)->set_min_refimenent_level(4);
       //tDatabase->get_background_mesh()->get_element(0)->put_on_refinement_queue();
       //tDatabase->flag_element( 0 );
       /*
       // refine the first element three times
       for( uint tLevel = 0; tLevel < 4; ++tLevel )
       {
           tDatabase->flag_element( 0 );

           // manually refine, do not reset pattern
          // tDatabase->get_background_mesh()->perform_refinement();
           tDatabase->perform_refinement( moris::hmr::RefinementMode::SIMPLE, false );
       }

       // update database etc
       tDatabase->perform_refinement( moris::hmr::RefinementMode::SIMPLE, false ); */

       //tDatabase->perform_refinement( moris::hmr::RefinementMode::LAGRANGE_REFINE, false );
       //tDatabase->perform_refinement( moris::hmr::RefinementMode::BSPLINE_REFINE, false );

       tHMR.flag_element( 0 );
       tHMR.get_database()->get_background_mesh()->get_element( 0 )->set_min_refimenent_level( 4 );

       tHMR.perform_refinement(  moris::hmr::RefinementMode::LAGRANGE_REFINE );
       tHMR.perform_refinement(  moris::hmr::RefinementMode::BSPLINE_REFINE );

       tHMR.finalize();
       moris::uint tBplineOrder = 2;
       moris::uint tLagrangeOrder = 2;
       auto tMesh = tHMR.create_mesh( tLagrangeOrder );
       uint tNumCoeffs = tMesh->get_num_coeffs( tBplineOrder );

       for( uint k=0; k<tNumCoeffs; ++k )
       {

       	std::string tLabel = "BSPline_" + std::to_string( k+1 );

       	auto tField = tMesh->create_field( tLabel, tBplineOrder );

   		Matrix<DDRMat> & tCoeffs = tField->get_coefficients();

   		tCoeffs.set_size( tMesh->get_num_coeffs( tBplineOrder ), 1, 0.0 );

   		tCoeffs( k ) = 1.0;

   		tField->evaluate_node_values();
       }



       //tHMR.perform_refinement_and_map_fields();

       tHMR.save_to_exodus( "Mesh1.exo" );
       tHMR.save_bsplines_to_vtk("Basis.vtk");

       //tHMR.save_last_step_to_exodus( "LastStep.exo" );

       //tHMR.save_to_hdf5( "Database.hdf5" );

       tHMR.save_coeffs_to_hdf5_file( "TMatrix.hdf5" );

//------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;

}

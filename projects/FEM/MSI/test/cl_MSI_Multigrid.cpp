
/*
 * cl_MSI_Multigrid.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "catch.hpp"

#include "typedefs.hpp"
#include "cl_Map.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "fn_equal_to.hpp"

#define protected public
#define private   public
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Node_Proxy.hpp"
#undef protected
#undef private

#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"

#include "cl_MTK_Mapper.hpp"

#include "fn_r2.hpp"

moris::real
LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

namespace moris
{
    namespace MSI
    {
    TEST_CASE("MSI_Multigrid","[MSI],[multigrid]")
    {
        if( moris::par_size() == 1 )
        {
             moris::hmr::Parameters tParameters;

             tParameters.set_number_of_elements_per_dimension( { { 2} , { 2 } } );
             tParameters.set_domain_offset( 0, 0 );

             uint tOrder = 2;
             tParameters.set_mesh_orders_simple( tOrder );
             tParameters.set_verbose( true );
             tParameters.set_multigrid( true );

             // create HMR object
             moris::hmr::HMR tHMR( tParameters );

             // flag first element
             tHMR.flag_element( 0 );
             tHMR.perform_refinement( moris::hmr::gRefinementModeBSpline );

             // finish mesh
             tHMR.finalize();

             // grab pointer to output field
             std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tOrder );

             // create field
             std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tOrder );
             std::shared_ptr< moris::hmr::Field > tExact = tMesh->create_field( "Exact", tOrder );

             // evaluate node values
             tField->evaluate_scalar_function( LevelSetFunction );
             tExact->get_node_values() = tField->get_node_values();

             // create mapper
             moris::mapper::Mapper tMapper( tMesh );

             // call mapping function
             tMapper.perform_mapping(
                     tField->get_label(),
                     EntityRank::NODE,
                     tField->get_label(),
                     tField->get_bspline_rank() );

             tField->evaluate_node_values();

             // save field to hdf5
             tField->save_field_to_hdf5("Circle.hdf5");

             // determine coefficient of determination
             moris::real tR2 = moris::r2( tExact->get_node_values(),
                     tField->get_node_values() );

             std::cout << "R2 " << tR2 << std::endl;
        }
    }
    }
}



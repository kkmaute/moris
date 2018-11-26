/*
 * cl_Dist_Vector_Test.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#ifdef MORIS_HAVE_PARALLEL
 #include "Epetra_MpiComm.h"
 #include <mpi.h>
#endif

#include "catch.hpp"

#include "fn_equal_to.hpp" // ALG/src

#include "typedefs.hpp" // COR/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Matrix_Vector_Factory.hpp" // DLA/src
#include "cl_Solver_Interface_Proxy.hpp" // DLA/src
#include "cl_Vector.hpp" // DLA/src

#include "MSI_Adof_Order_Hack.hpp"

#define protected public
#define private   public
#include "cl_MSI_Multigrid.hpp"
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

#include "cl_FEM_Node_Base.hpp"
#include "cl_FEM_Element.hpp"

#include "cl_MTK_Mapper.hpp"
#include "cl_FEM_IWG_L2.hpp"

#include "fn_r2.hpp"

namespace moris
{
TEST_CASE("DLA_Multigrid","[DLA],[DLA_multigrid]")
{
    if( moris::par_size() == 1 )
    {
        // order for this example
        moris::uint tOrder = 1;
        moris::MSI::gAdofOrderHack = tOrder;

        // create parameter object
        moris::hmr::Parameters tParameters;
        tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );
        tParameters.set_verbose( false );
        tParameters.set_multigrid( true );
        tParameters.set_bspline_truncation( true );
        tParameters.set_mesh_orders_simple( tOrder );

        // create HMR object
        moris::hmr::HMR tHMR( tParameters );

        // flag first element for refinement
        tHMR.flag_element( 0 );
        tHMR.perform_refinement( moris::hmr::gRefinementModeBSpline );

        tHMR.flag_element( 0 );
        tHMR.perform_refinement( moris::hmr::gRefinementModeBSpline );

        tHMR.finalize();

         // grab pointer to output field
         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tOrder );

         //tHMR.save_bsplines_to_vtk("BSplines.vtk");

         moris::map< moris::moris_id, moris::moris_index > tMap;
         tMesh->get_adof_map( tOrder, tMap );
         //tMap.print("Adof Map");

         //-------------------------------------------------------------------------------------------

         // create IWG object
         fem::IWG_L2 * tIWG = new moris::fem::IWG_L2( );

         map< moris_id, moris_index >   tCoefficientsMap;
         Cell< fem::Node_Base* >        tNodes;
         Cell< MSI::Equation_Object* >  tElements;

         // get map from mesh
         tMesh->get_adof_map( tOrder, tCoefficientsMap );

         // ask mesh about number of nodes on proc
         luint tNumberOfNodes = tMesh->get_num_nodes();

         // create node objects
         tNodes.resize( tNumberOfNodes, nullptr );

         for( luint k = 0; k < tNumberOfNodes; ++k )
         {
             tNodes( k ) = new fem::Node( &tMesh->get_mtk_vertex( k ) );
         }

         // ask mesh about number of elements on proc
         luint tNumberOfElements = tMesh->get_num_elems();

         // create equation objects
         tElements.resize( tNumberOfElements, nullptr );

         for( luint k=0; k<tNumberOfElements; ++k )
         {
             // create the element
             tElements( k ) = new fem::Element( & tMesh->get_mtk_cell( k ),
                                                tIWG,
                                                tNodes );
         }

         MSI::Model_Solver_Interface * tMSI = new moris::MSI::Model_Solver_Interface( tElements,
                                                                                      tMesh->get_communication_table(),
                                                                                      tCoefficientsMap,
                                                                                      tMesh->get_num_coeffs( tOrder ),
                                                                                      tMesh.get() );

         moris::Matrix< DDSMat > tExternalIndices( 9, 1 );
         tExternalIndices( 0, 0 ) = 17;
         tExternalIndices( 1, 0 ) = 18;
         tExternalIndices( 2, 0 ) = 21;
         tExternalIndices( 3, 0 ) = 22;
         tExternalIndices( 4, 0 ) = 23;
         tExternalIndices( 5, 0 ) = 25;
         tExternalIndices( 6, 0 ) = 27;
         tExternalIndices( 7, 0 ) = 26;
         tExternalIndices( 8, 0 ) = 28;

         moris::Matrix< DDSMat > tInternalIndices;

         tMSI->read_multigrid_maps( 2, tExternalIndices, 0, tInternalIndices );

         CHECK( equal_to( tInternalIndices( 0, 0 ), 0 ) );
         CHECK( equal_to( tInternalIndices( 1, 0 ), 1 ) );
         CHECK( equal_to( tInternalIndices( 2, 0 ), 2 ) );
         CHECK( equal_to( tInternalIndices( 3, 0 ), 3 ) );
         CHECK( equal_to( tInternalIndices( 4, 0 ), 4 ) );
         CHECK( equal_to( tInternalIndices( 5, 0 ), 5 ) );
         CHECK( equal_to( tInternalIndices( 6, 0 ), 6 ) );
         CHECK( equal_to( tInternalIndices( 7, 0 ), 7 ) );
         CHECK( equal_to( tInternalIndices( 8, 0 ), 8 ) );

         delete tMSI;
         delete tIWG;

         for( luint k=0; k<tNumberOfElements; ++k )
         {
             // create the element
             delete tElements( k );
         }

         for( luint k = 0; k < tNumberOfNodes; ++k )
         {
             delete tNodes( k );
         }

    }
}

}



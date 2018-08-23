/*
 * main.cpp
 *
 *  Created on: Dec 22, 2017
 *      Author: doble
 */

// MORIS header files.
#include <fn_r2.hpp>
#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_HMR.hpp" // HMR/src
#include "cl_FEM_Element.hpp" // FEM/INT/src


#include "cl_Model_Solver_Interface.hpp"
#include "cl_Equation_Object.hpp"
#include "cl_FEM_IWG_L2_Test.hpp"
#include "cl_MSI_Node.hpp"

#include "cl_MDL_Model.hpp"





moris::Comm_Manager gMorisComm;

using namespace moris;

real
distance( const Mat< real > & aPoint )
{
    real tValue = 0;
    uint tNumberOfDimensions = aPoint.length();

    for( uint k=0; k< tNumberOfDimensions; ++k )
    {
        tValue += std::pow( aPoint( k ), 2 );
    }
    return std::sqrt( tValue );
}

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager(&argc, &argv);

    // print welcome banner and system information
    moris::print_banner( argc, argv );

//------------------------------------------------------------------------------
    // In this part, an HMR object is created

    // create settings object
    /* ( this is actually my internal settings object. As far as I understood,
                 the user is supposed to control HMR using a parameter list.
                 This is not implemented yet ) */

    hmr::Parameters tParameters;

    // create a Mat for a 2D object
    //Mat< luint > tNumberOfElements = { { 2 }, { 2 }, { 2 } };
    Mat< luint > tNumberOfElements = { { 2 }, { 2 } };
    // pass number of elements to settings
    tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

    // make mesh output silent
    tParameters.set_verbose( true );

    // set maximum interpolation degree
    tParameters.set_mesh_orders_simple( 1 );

    // buffer size must be set at least to max polynomial if truncation is used
    //tParameters.set_buffer_size( tParameters.get_max_polynomial() );

    tParameters.set_buffer_size( 0 );
    tParameters.set_bspline_truncation( false );

    hmr::HMR tHMR( &tParameters );

    // create a mesh interface
    auto  tMesh = tHMR.create_interface();

    // < * usually, this is where the refinement logic would happen * >

//------------------------------------------------------------------------------

    // create IWG object
    moris::fem::IWG_L2_Test tIWG;

    auto tField = tHMR.create_field( "Field", 0 );

    //for( uint tLevel = 1; tLevel<4; ++tLevel )
    {

        uint tNumberOfElements = tHMR.get_number_of_elements_on_proc();

        for( uint e=0; e<tNumberOfElements; e = e+3 )
        {
            tHMR.flag_element( e );
        }
        tHMR.perform_refinement();

        uint tOrder = 1;
        //for( uint tOrder = 1; tOrder <= 3; ++tOrder )
        {
            auto tBlock = tMesh.get_block_by_index( tOrder-1 );


            // how many cells exist on current proc
            auto tNumberOfCells = tBlock->get_number_of_cells();

            for( luint k=0; k<tNumberOfCells; ++k )
            {
                tBlock->get_cell_by_index( k )->set_t_matrix_flag();
            }

            // after all equation objects  are created, calculate the T-Matrices
            tMesh.finalize();


            auto tNumberOfNodes = tBlock->get_number_of_vertices();

            // initialize cell
            moris::Cell< moris::MSI::Equation_Object* > tListEqnObj( tNumberOfCells, nullptr );

            Mat< real > tWeakBCs( tNumberOfNodes, 1 );
            for( uint k=0;  k<tNumberOfNodes; ++k )
            {
                tWeakBCs( k ) = distance( tBlock->get_vertex_by_index( k )->get_coords() );
            }


            Mat< real > tResult;
            mdl::Model tModel( tMesh, tIWG, tWeakBCs, tResult );


            // create matrix with node values
            //Mat< real > tNodeValues( tBlock->get_number_of_vertices(), 1, 0.0 );
            Mat< real > & tNodeValues( tField->get_data() );
            tNodeValues.set_size( tBlock->get_number_of_vertices(), 1, 0.0 );
            // copy node values from equation object
            for ( auto tElement : tModel.mEquationObjects )
            {
                tElement->get_pdof_values( tNodeValues );
            }

            tNodeValues.print("Values");

            //tHMR.add_field( "Field", tOrder, tNodeValues );


//------------------------------------------------------------------------------


            // compare output with exact solution
            // https://en.wikipedia.org/wiki/Coefficient_of_determination

            Mat< real > tVertexNorm( tNumberOfNodes, 1 );


            //real tAvg = 0;

            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                // get pointer to vertex
                auto tVertex = tBlock->get_vertex_by_index( k );

                tVertexNorm( k ) = norm ( tVertex->get_coords() );
            }


            real tR2 = r2( tNodeValues, tVertexNorm );

            std::cout << "Coefficient of determination: " << tR2<< std::endl;

        }
    }

    // get number of vertices
    tHMR.save_to_exodus( "Mesh.exo" );

//------------------------------------------------------------------------------
    // delete iwg pointer
    //delete tIWG;

//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

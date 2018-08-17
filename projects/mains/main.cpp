/*
 * main.cpp
 *
 *  Created on: Dec 22, 2017
 *      Author: doble
 */

// MORIS header files.
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
moris::Comm_Manager gMorisComm;

using namespace moris;

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
    Mat< luint > tNumberOfElements = { { 2 }, { 2 }, { 2 }  };

    // pass number of elements to settings
    tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

    // make mesh output silent
    tParameters.set_verbose( false);

    // set maximum interpolation degree
    tParameters.set_max_polynomial( 2 );

    // buffer size must be set to max polynomial if truncation is used
    tParameters.set_buffer_size( tParameters.get_max_polynomial() );
    tParameters.set_bspline_truncation( true );

    hmr::HMR tHMR( &tParameters );

    // < * usually, this is where the refinement logic would happen * >
    tHMR.flag_element( 0 );
    tHMR.flag_element( 7 );
    tHMR.perform_refinement();

//------------------------------------------------------------------------------

    // create a mesh interface
    auto tMesh = tHMR.create_interface();

    auto tBlock = tMesh.get_block_by_index( 1 );


    // create pointer to IWG object
    auto tIWG = new moris::fem::IWG_L2_Test();

    //auto tElement = new fem::Element( tBlock->get_cell_by_index( 0 ), tIWG );

    // how many cells exist on current proc
    auto tNumberOfCells = tBlock->get_number_of_cells();

    // initialize cell
    moris::Cell< moris::MSI::Equation_Object* > tListEqnObj( tNumberOfCells, nullptr );

    std::cout << "Number of Cells: " << tNumberOfCells << std::endl;

    // populate cell
    for( luint k=0; k<tNumberOfCells; ++k )
    {
        tListEqnObj( k ) = new moris::MSI::Equation_Object(
                tBlock->get_cell_by_index( k ),
                tIWG );
    }

    // after all equation objects  are created, calculate the T-Matrices
    tMesh.finalize();

    // return the communication table
    moris::Mat< moris::uint > tCommTable = tMesh.get_communication_table();

    tCommTable.print("CommunucationTable");

//------------------------------------------------------------------------------
    if( par_size() == 1)
    {
        // this part does not work yet in parallel
        moris::MSI::Model_Solver_Interface tMSI( tListEqnObj, tCommTable );

        tMSI.solve_system( tListEqnObj );
    }
//------------------------------------------------------------------------------

    // create matrix with node values
    Mat< real > tNodeValues( tBlock->get_number_of_vertices(), 1 );

    // copy node values from equation object
    for ( auto tElement : tListEqnObj )
    {
        tElement->get_pdof_values( tNodeValues );
    }

    //tNodeValues.print("Values");

    tHMR.add_field( "Field", tOrder, tNodeValues );

    // clean up memory
    for ( auto tElement : tListEqnObj )
    {
        tElement->get_pdofs_values();
    }

    // clean up memory
    for ( auto tElement : tListEqnObj )
    {
        delete tElement;
    }
    tHMR.save_to_exodus( tOrder, "Mesh.exo");


//------------------------------------------------------------------------------
    // delete iwg pointer
    delete tIWG;

//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

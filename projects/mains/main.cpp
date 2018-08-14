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
    Mat< luint > tNumberOfElements = { { 1 }, { 1 } };

    // pass number of elements to settings
    tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

    // make mesh output silent
    tParameters.set_verbose( false );

    // set maximum interpolation degree
    tParameters.set_max_polynomial( 3 );

    // buffer size must be set to max polynomial if truncation is used
    tParameters.set_buffer_size( tParameters.get_max_polynomial() );
    tParameters.set_bspline_truncation( true );

    /* initialize hmr object
               using the present settings, HMR will create three B-Spline meshes
               and three Lagrange meshes. There will be an option that allows
               selecting the mesh orders individually */
    hmr::HMR tHMR( &tParameters );

    // < * usually, this is where the refinement logic would happen * >

//------------------------------------------------------------------------------
    // This is to be placed in the model

    // create a mesh interface
    auto tMesh = tHMR.create_interface();

    // print number of meshes
    // std::cout << " number of blocks " << tMesh.get_number_of_blocks() << std::endl;

    // pick second order block on mesh
    auto tBlock = tMesh.get_block_by_index( 1 );

    // create element out of first element on mesh
/*    auto tElement = new fem::Element( tBlock->get_cell_by_index( 0 ) );

    // after all elements are created, calculate the T-Matrices
    tMesh.finalize();

    // this tells how many nodes exists on the block
    //auto tNumberOfNodes = tBlock->get_number_of_vertices();

    // get vertex of block
    // auto tVertex = tBlock->get_vertex_by_index( 0 );

    // this tells how many elements exist in the block
    //auto tNumberOfElements = tBlock->get_number_of_cells();

    // calc mass matrix for element
    Mat< real > tM;
    tElement->eval_mass( tM );

    tM.print("Mass");


    // how many nodes belong to this element
    auto tNumberOfNodes = tElement->get_number_of_nodes();

    Cell< mtk::Vertex* > tNodes = tElement->get_vertex_pointers();

    for( uint k=0; k<tNumberOfNodes; ++k )
    {
        // get pointers to DOFs that are connected to this element
        // auto tDOFs = tNodes( k )->get_adof_pointers();

        // get ids of B-Splines that are connected to this element
        auto tIDs = tNodes( k )->get_adof_ids();

        // get T-Matrix of this node
        auto tTMatrix = tNodes( k )->get_t_matrix();

        // get ownners of the dofs
        //auto tDofOwners = tNodes( k )->get_adof_owners();

        std::cout << "Node " << k << std::endl;

        tIDs.print("DOF IDs");

        tTMatrix->print("T-Matrix");

        //tDofOwners.print("DOF Owner");

    } */

	//get number of elements from mesh
	auto tNumberOfCells = tBlock->get_number_of_cells();

	// initialize cell
    moris::Cell< moris::MSI::Equation_Object* > tListEqnObj( tNumberOfCells, nullptr );

    // populate cell
    for( luint k=0; k<tNumberOfCells; ++k )
    {
       tListEqnObj( k ) = new moris::MSI::Equation_Object( tBlock->get_cell_by_index( k ) );
    }

    // after all equation objects  are created, calculate the T-Matrices
    tMesh.finalize();


//-------------------------------------------------------------------
    moris::uint tNumEquationObjects = tListEqnObj.size();
    moris::MSI::Model_Solver_Interface tMSI( tNumEquationObjects, tListEqnObj );

    tMSI.solve_system();

//------------------------------------------------------------------------------

    // clean up memory
    for ( auto tElement : tListEqnObj )
    {
    	delete tElement;
    }


    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

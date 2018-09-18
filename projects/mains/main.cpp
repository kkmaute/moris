/*
 * main.cpp
 *
 *  Created on: Dec 22, 2017
 *      Author: doble
 */

// MORIS header files.

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src

#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Refinement_Manager.hpp"

#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

moris::Comm_Manager gMorisComm;

using namespace moris;
using namespace hmr;

// a simple function for the levelset
real
LevelSet( const Mat< real > & aPoint )
{
    return norm( aPoint ) - 1.5;
}


int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // print welcome banner and system information
    //moris::print_banner( argc, argv );

    /*
     * Step 1: Load 2D mesh From file
     */
    HMR tHMR( "Tutorial1.hdf5" );

    // fixme: after initialization, remember to calculate T-Matrices
    // for input meshes

    //------------------------------------------------------------------------------
    /*
     * Step 2: Create an MTK Field and load values from binary file
     */

    /*
     * create an MTK Mesh pointer for the input pattern
     */
    Mesh * tInputMesh = tHMR.create_input_mesh();

    /*
     * create a field object that contains mesh data
     */
    mtk::Field * tInputField = tInputMesh->create_field( "LevelSet" );

    /*
     * load node values from binary file
     */
    load_matrix_from_binary_file( * tInputField->get_node_values(), "Tutorial1.bin" );

    /* Step 3: Create an MTK Refinement Manager and Identify
     *         Volume and Surface Cells, and flag them for refinement
     */

    /*
     * Matrix with indices of cells where all level set values are < 0
     */
    moris::Mat< moris::moris_index > tVolumeCells;

    /*
     * Matrix with indices of intersected cells
     */
    moris::Mat< moris::moris_index > tSurfaceCells;

    /*
     * Create a refinement Manager Object
     */
    moris::mtk::Refinement_Manager tRefMan;

    // perform refinement check
    tRefMan.find_volume_and_surface_cells(
            tVolumeCells,
            tSurfaceCells,
            tInputField );

    tSurfaceCells.print("Surface");
    tVolumeCells.print("Volume");

    delete tInputField;
    delete tInputMesh;

    tHMR.save_to_exodus( "Mesh.exo");

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

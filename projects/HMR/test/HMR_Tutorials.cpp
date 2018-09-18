/*
 * HMR_Tutorials.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: messe
 */

#include "catch.hpp"





#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Refinement_Manager.hpp"

#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

using namespace moris;
using namespace hmr;

// a simple function for the levelset
real
LevelSet( const Mat< real > & aPoint )
{
    return norm( aPoint ) - 1.5;
}


TEST_CASE("HMR Tutorial","[HMR_Tutorial]")
{
//------------------------------------------------------------------------------
    if( moris::par_rank() == 0 )
    {
        /*
         * This example creates a simple 2x2 mesh in serial and refines it according
         * to a given LevelSet
         */
        SECTION( "Initialize Simple 2D Mesh" )
        {

//------------------------------------------------------------------------------
            /*
             * Step 1: create a parameter list
             */

            /*
             * Create a default parameter list. It already as default settings
             * which are not of interest at the moment.
             */
            ParameterList tParameters = create_parameter_list();

            /*
             * the coarsest grid will have 3 elements in x-direction,
             * and 2 elements in y-direction
             */
            tParameters.set("number_of_elements_per_dimension", "3, 2" );


            /*
             * the domain will have a width of 3 and a height of 2 units
             */
            tParameters.set("domain_dimensions", "3, 2" );

            /*
             * set the coordinates of the lower left node on the mesh
             */
            tParameters.set("domain_offset", "0, 0");

            /*
             * create a second order mesh
             */
            tParameters.set( "interpolation_order", "2" );

//------------------------------------------------------------------------------

            /*
             * Step 2: HMR object
             */

            /*
             * create the HMR main object from the parameters
             */
            HMR tHMR( tParameters );


//------------------------------------------------------------------------------

            /*
             * Step 3: Create an MTK Field and fill it with nodal values
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
             * evaluate levelset on input mesh
             */
            tInputField->evaluate_scalar_function( LevelSet );

//------------------------------------------------------------------------------

            /*
             * Step 4: Create an MTK Refinement Manager and Identify
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


            /*
             * Flag volume elements for refinement.
             * The Refinement_Manager does not know the concept of levels.
             * A level lowpass can be specified as follows
             */
            tHMR.flag_elements( tVolumeCells,
                    tParameters.get< moris::sint >( "max_volume_refinement_level") );

            /*
             * Flag surface elements for refinement.
             * The Refinement_Manager does not know the concept of levels.
             * A level lowpass can be specified as follows
             */
            tHMR.flag_elements( tSurfaceCells,
                    tParameters.get< moris::sint >( "max_surface_refinement_level") );

            /*
             * with all elements flagged, we can now tell HMR to perform the refinement routine
             */
            tHMR.perform_refinement();

//------------------------------------------------------------------------------

            /*
             * Step 5: Map field to output mesh
             */

            /*
             * before the mapping can be performed, the input and the output patterns
             * must be united
             */
            tHMR.create_union_pattern();

            /*
             * create a pointer to the output mesh
             */
            Mesh * tOutputMesh = tHMR.create_output_mesh();

            /*
             * create an output field and copy label from input field
             */
            mtk::Field * tOutputField = tHMR.map_field_on_mesh( tInputField, tOutputMesh );

//------------------------------------------------------------------------------

            /*
             * Step 6: Write data
             */

            /*
             * Write mesh to file
             */
            tHMR.save_to_hdf5( "Tutorial_1.hdf5" );

            /*
             * Write Output field to file
             */
            save_matrix_to_binary_file( *tOutputField->get_node_values(), "Tutorial_1.bin" );

            /*
             * STK_Implementation is currently not linked to the new MTK interface.
             * This is a workaround for now.
             */
            tHMR.save_to_exodus("Tutorial_1.exo");

//------------------------------------------------------------------------------

            /*
             * Step 7: Tidy up memory
             */

            /*
             * delete pointers: Output and input fields and meshes must be deleted
             */
            delete tInputField;
            delete tInputMesh;

            delete tOutputField;
            delete tOutputMesh;

//------------------------------------------------------------------------------
        }

        SECTION( "Load 2D Mesh from File and refine again" )
        {
            /*
             * Step 1: Load 2D mesh From file
             */
            HMR tHMR( "Tutorial_1.hdf5" );

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
            load_matrix_from_binary_file( * tInputField->get_node_values(), "Tutorial_1.bin" );


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
            /*
             * Flag volume elements for refinement.
             * The Refinement_Manager does not know the concept of levels.
             * A level lowpass can be specified as follows
             */
            tHMR.flag_elements( tVolumeCells,
                    tHMR.get_parameters()->get_max_volume_level() );

            /*
             * Flag surface elements for refinement.
             * The Refinement_Manager does not know the concept of levels.
             * A level lowpass can be specified as follows
             */
            tHMR.flag_elements( tSurfaceCells,
                   tHMR.get_parameters()->get_max_surface_level() );

            /*
             * with all elements flagged, we can now tell HMR to perform the refinement routine
             */
            tHMR.perform_refinement();

//------------------------------------------------------------------------------

            /*
             * Step 4: Map field to output mesh
             */

            /*
             * before the mapping can be performed, the input and the output patterns
             * must be united
             */
            tHMR.create_union_pattern();

            /*
             * create a pointer to the output mesh
             */
            Mesh * tOutputMesh = tHMR.create_output_mesh();

            /*
             * create an output field and copy label from input field
             */
            mtk::Field * tOutputField = tHMR.map_field_on_mesh( tInputField, tOutputMesh );

//------------------------------------------------------------------------------

            /*
             * Step 5: Write data
             */

            /*
             * Write mesh to file
             */
            tHMR.save_to_hdf5( "Tutorial_2.hdf5" );

            /*
             * Write Output field to file
             */
            save_matrix_to_binary_file( *tInputField->get_node_values(), "Tutorial_2.bin" );

            /*
             * STK_Implementation is currently not linked to the new MTK interface.
             * This is a workaround for now.
             */
            tHMR.save_to_exodus("Tutorial_2.exo");

//------------------------------------------------------------------------------

            /*
             * Step 6: Tidy up memory
             */

            /*
             * delete pointers: Output and input fields and meshes must be deleted
             */
            delete tInputField;
            delete tInputMesh;

            delete tOutputField;
            delete tOutputMesh;
        }
    }

}

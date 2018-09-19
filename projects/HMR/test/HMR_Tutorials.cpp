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
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

using namespace moris;
using namespace hmr;

/*!
 * the folowing function is used as level set for the tutorial
 *
 * \code{.cpp}
 * real
 * LevelSetFunction( const Mat< real > & aPoint )
 * {
 *     return norm( aPoint ) - 1.2;
 * }
 * \endcode
 */
real
LevelSetFunction( const Mat< real > & aPoint )
{
    return norm( aPoint ) - 1.2;
}


TEST_CASE("HMR Tutorial","[HMR_Tutorial]")
{
//------------------------------------------------------------------------------
    if( moris::par_rank() == 0 )
    {
        /*!
         * <h1>Tutorial 1: Initialize a Simple 2D Mesh"</h1>
         * This example creates a simple 2x3 mesh in serial and refines it according
         * to a given LevelSet.
         */
        SECTION( "Initialize a Simple 2D Mesh" )
        {

//------------------------------------------------------------------------------
            /*!
             * <b> Step 1: create a parameter list </b>
             */

            /*!
             * Create a default parameter list. It already as default settings
             * which are not of interest at the moment.
             *
             * \code{.cpp}
             * ParameterList tParameters = create_parameter_list();
             * \endcode
             */
            ParameterList tParameters = create_parameter_list();

            /*!
             * the coarsest grid will have 3 elements in x-direction,
             * and 2 elements in y-direction
             *
             * \code{.cpp}
             * tParameters.set("number_of_elements_per_dimension", "3, 2" );
             * \endcode
             */
            tParameters.set("number_of_elements_per_dimension", "3, 2" );


            /*!
             * the domain will have a width of 3 and a height of 2 units
             *
             * \code{.cpp}
             * tParameters.set("domain_dimensions", "3, 2" );
             * \endcode
             */
            tParameters.set("domain_dimensions", "3, 2" );

            /*!
             * set the coordinates of the lower left node on the mesh
             *
             * \code{.cpp}
             * tParameters.set("domain_offset", "0, 0");
             * \endcode
             */
            tParameters.set("domain_offset", "0, 0");

            /*!
             * Create a second order mesh. The parameter has to be set as
             * as string. In the future, HMR will also support
             * multiple interpolation orders, such as "2, 3".
             *
             * \code{.cpp}
             * tParameters.set( "interpolation_order", "2" );
             * \endcode
             */
            tParameters.set( "interpolation_order", "2" );

//------------------------------------------------------------------------------

            /*!
             * <b> Step 2: HMR object </b>
             */

            /*!
             * create the HMR main object from the parameters
             *
             * \code{.cpp}
             * HMR tHMR( tParameters );
             * \endcode
             */
            HMR tHMR( tParameters );


//------------------------------------------------------------------------------

            /*!
             * <b> Step 3: Create an MTK Field and fill it with nodal values </b>
             */

            /*!
             * The hmr::Mesh pointer provides an interface to the MTK API.
             * A pointer needs to be created for both the input and the
             * output mesh.
             *
             * \code{.cpp}
             *  Mesh * tInputMesh = tHMR.create_input_mesh();
             * \endcode
             */
            Mesh * tInputMesh = tHMR.create_input_mesh();


            /*!
             * Create a field object that is linked to the input mesh.
             * Name it "LevelSet".
             *
             * \code{.cpp}
             * mtk::Field * tInputField = tInputMesh->create_field( "LevelSet" );
             * \endcode
             */
            mtk::Field * tInputField = tInputMesh->create_field( "LevelSet" );

            /*!
             * The function LevelSetFunction, is evaluated on the nodes.
             * Note that this call does not initialize B-Spline coefficients.
             *
             * \code{.cpp}
             * tInputField->evaluate_scalar_function( LevelSetFunction );
             * \endcode
             */
            tInputField->evaluate_scalar_function( LevelSetFunction );

//------------------------------------------------------------------------------

            /*
             * <b> Step 4: Create an MTK Refinement Manager and Identify
             *         Volume and Surface Cells, and flag them for refinement </b>
             */


            /*!
             * the following lines do the following
             * - create a geometry engine
             * - ask HMR which elements can be considered for the volume refinement
             * - ask the geometry engine for cells enclosed in the level set
             * - tell HMR to put these elements on the waiting list for refinement
             * - ask HMR which elements can be considered for the surface refinement
             * - ask the geometry engine for cells that are intersected by the level set
             * - tell HMR to put these elements on the waiting list for refinement
             *
             * These steps can be summarized by calling
             *
             * \code{.cpp}
             * tHMR.flag_volume_and_surface_elements( tInputField );
             * \endcode
             *
             * The individual steps of this funciton are as follows:
             */


            /*!
             * create a geomery engine
             *
             * \code{.cpp}
             * gen::Geometry_Engine tRefinementManager;
             * \endcode
             */
            gen::Geometry_Engine tRefinementManager;

            /*!
             * create a Cell of mtk::Cells that are to be considered
             * \code{.cpp}
             * Cell< mtk::Cell* > tCandidates;
             * \endcode
             */
            Cell< mtk::Cell* > tCandidates;

            /*!
             * create a Cell of mtk::Cells that are identified by the Geometry Engine
             *
             * \code{.cpp}
             * Cell< mtk::Cell* > tRefinementList;
             * \endcode
             */
            Cell< mtk::Cell* > tRefinementList;

            /*!
             * ask HMR for candidates for the volume refinement
             *
             * \code{.cpp}
             * tHMR.get_candidates_for_volume_refinement( tCandidates );
             * \endcode
             */
            tHMR.get_candidates_for_volume_refinement( tCandidates );

            /*!
             * investigate these elements by the Geometry Engine
             *
             * \code{.cpp}
             * tRefinementManager.find_cells_within_levelset(
             *       tRefinementList,
             *       tCandidates,
             *       aScalarField );
             * \endcode
             */
            tRefinementManager.find_cells_within_levelset(
                    tRefinementList,
                    tCandidates,
                    tInputField );

            /*!
             * add these elements to the waiting list in HMR
             *
             *
             * \code{.cpp}
             * tHMR.flag_elements( tRefinementList );
             * \endcode
             */
            tHMR.flag_elements( tRefinementList );

            /*!
             * ask HMR for candidates for the surface refinement
             *
             * \code{.cpp}
             * tHMR.get_candidates_for_surface_refinement( tCandidates );
             * \endcode
             */
            tHMR.get_candidates_for_surface_refinement( tCandidates );


            /*!
             * investigate these elements by the Geometry Engine
             *
             * \code{.cpp}
             * tRefinementManager.find_cells_intersected_by_levelset(
             *       tRefinementList,
             *       tCandidates,
             *       aScalarField );
             * \endcode
             */
            tRefinementManager.find_cells_intersected_by_levelset(
                    tRefinementList,
                    tCandidates,
                    tInputField );

            /*!
             * add these elements to the waiting list in HMR
             *
             *
             * \code{.cpp}
             * tHMR.flag_elements( tRefinementList );
             * \endcode
             */
            tHMR.flag_elements( tRefinementList );

            /*!
             * The process of flagging elements can be repeated for an arbitrary number
             * of times. Once all elements are flagged, the refinement procedure is
             * called.
             *
             * \code{.cpp}
             * tHMR.perform_refinement();
             * \endcode
             */
            tHMR.perform_refinement();

//------------------------------------------------------------------------------

            /*!
             * <b> Step 5: Map field to output mesh </b>
             *
             * Now we have a mesh that is refined to the specified criteria.
             * The input field must now be projected onto the output.
             */

            /*!
             * create a pointer to the output mesh
             *
             * \code{.cpp}
             * Mesh * tOutputMesh = tHMR.create_output_mesh();
             * \endcode
             */
            Mesh * tOutputMesh = tHMR.create_output_mesh();

            /*!
             * Call the L2 projection and write pointer to output mesh
             *
             * \code{.cpp}
             * mtk::Field * tOutputField = tHMR.map_field_on_mesh( tInputField, tOutputMesh );
             * \endcode
             */
            mtk::Field * tOutputField = tHMR.map_field_on_mesh( tInputField, tOutputMesh );

//------------------------------------------------------------------------------

            /*!
             * <b> Step 6: Write data </b>
             *
             * Now we can save the mesh into a file.
             */

            /*!
             * The HDF5 file stores the settings and the refinement pattern.
             *
             * \code{.cpp}
             * tHMR.save_to_hdf5( "Tutorial_1.hdf5" );
             * \endcode
             */
            tHMR.save_to_hdf5( "Tutorial_1.hdf5" );

            /*!
             * In order to save the field, we grab the pointer to node values of the mesh
             * and save it to a binary file. Alternatively, one could also store the
             * B-Spline coefficients the same way.
             *
             * \code{.cpp}
             * save_matrix_to_binary_file( *tOutputField->get_node_values(), "Tutorial_1.bin" );
             * \endcode
             */
            save_matrix_to_binary_file( *tOutputField->get_node_values(), "Tutorial_1.bin" );

            /*!
             * STK_Implementation is currently not fully linked to the new MTK interface.
             * HMR's functionality to store exodus files is a temporary solution
             * until this issue is resolved.
             *
             * \code{.cpp}
             * tHMR.save_to_exodus("Tutorial_1.exo");
             * \endcode
             */
            tHMR.save_to_exodus("Tutorial_1.exo");

//------------------------------------------------------------------------------

            /*!
             * <b> Step 7: Tidy up memory </b>
             *
             * Finally, we can delete the pointers for the meshes and the fields.
             */

            /*
             * delete pointers: Output and input fields and meshes must be deleted
             * \code{.cpp}
             * delete tInputField;
             * delete tInputMesh;
             * delete tOutputField;
             * delete tOutputMesh;
             * \endcode
             */
            delete tInputField;
            delete tInputMesh;

            delete tOutputField;
            delete tOutputMesh;

//------------------------------------------------------------------------------
        }

        /*!
         * <h1>Tutorial 2 : Load 2D Mesh from File and refine again"</h1>
         * In the following example, a mesh and a nodal field are loaded
         * from existing files. The refinement is performed and
         * the new mesh is storen into an output.
         */
        SECTION( "Load 2D Mesh from File and refine again" )
        {
            /*!
             * <b> Step 1: Load 2D mesh From file </b>
             *
             * \code{.cpp}
             * HMR tHMR( "Tutorial_1.hdf5" );
             * \endcode
             */
            HMR tHMR( "Tutorial_1.hdf5" );

//------------------------------------------------------------------------------
            /*!
             * <b> Step 2: Create an MTK Field and load values from binary file </b>
             */

            /*!
             * create an MTK Mesh pointer for the input pattern
             * \code{.cpp}
             * Mesh * tInputMesh = tHMR.create_input_mesh();
             * \endcode
             */
            Mesh * tInputMesh = tHMR.create_input_mesh();

            /*!
             * create a field object that contains mesh data
             *
             * \code{.cpp}
             * mtk::Field * tInputField = tInputMesh->create_field( "LevelSet" );
             * \endcode
             */
            mtk::Field * tInputField = tInputMesh->create_field( "LevelSet" );

            /*!
             * load node values from binary file
             *
             * \code{.cpp}
             * load_matrix_from_binary_file( * tInputField->get_node_values(), "Tutorial_1.bin" );
             * \endcode
             */

            // fixme: Make input mesh calculate T-Matrices after initialization.
            //        Then, coefficents can be loaded from binary as well
            load_matrix_from_binary_file( * tInputField->get_node_values(), "Tutorial_1.bin" );


            /*! <b> Step 3: Refine the mesh </b>
             *
             * The task of creating a refinement manager and flagging elements
             * according a surface and volume criterion is used frequently by HMR.
             * Therefore, the steps from step 4 of the first example can be summarized
             * by the following lines
             *
             * \code{.cpp}
             * tHMR.flag_volume_and_surface_elements( tInputField );
             * tHMR.perform_refinement();
             * \endcode
             */

            // flag volume and surface elements
            tHMR.flag_volume_and_surface_elements( tInputField );

            // refine mesh
            tHMR.perform_refinement();

//------------------------------------------------------------------------------

            /*!
             * <b> Step 4: Map field to output mesh </b>
             */

            /*!
             * create a pointer to the output mesh
             * \code{.cpp}
             * Mesh * tOutputMesh = tHMR.create_output_mesh();
             * \endcode
             */
            Mesh * tOutputMesh = tHMR.create_output_mesh();

            /*!
             * create an output field and copy label from input field
             * \code{.cpp}
             * mtk::Field * tOutputField = tHMR.map_field_on_mesh( tInputField, tOutputMesh );
             * \endcode
             */
            mtk::Field * tOutputField = tHMR.map_field_on_mesh( tInputField, tOutputMesh );

//------------------------------------------------------------------------------

            /*!
             * <b> Step 5: Write data </b>
             */

            /*!
             * Write mesh to file
             * \code{.cpp}
             * tHMR.save_to_hdf5( "Tutorial_2.hdf5" );
             * \endcode
             */
            tHMR.save_to_hdf5( "Tutorial_2.hdf5" );

            /*!
             * Write output field to file
             * \code{.cpp}
             * save_matrix_to_binary_file( *tInputField->get_node_values(), "Tutorial_2.bin" );
             * \endcode
             */
            save_matrix_to_binary_file( *tInputField->get_node_values(), "Tutorial_2.bin" );

            /*!
             * STK_Implementation is currently not linked to the new MTK interface.
             * This is a workaround for now.
             * \code{.cpp}
             * tHMR.save_to_exodus("Tutorial_2.exo");
             * \endcode
             */
            tHMR.save_to_exodus("Tutorial_2.exo");

//------------------------------------------------------------------------------

            /*!
             * <b> Step 6: Tidy up memory </b>
             */

            /*!
             * delete pointers: Output and input fields and meshes must be deleted
             * \code{.cpp}
             * delete tInputField;
             * delete tInputMesh;
             * delete tOutputField;
             * delete tOutputMesh;
             * \endcode
             */
            delete tInputField;
            delete tInputMesh;

            delete tOutputField;
            delete tOutputMesh;
        }
    }

}

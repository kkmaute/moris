//------------------------------------------------------------------------------

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_MTK_Field.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------
// HMR
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

//------------------------------------------------------------------------------

// select namespaces
using namespace moris;
using namespace hmr;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;

//------------------------------------------------------------------------------

/*!
 * the folowing function is used as level set for the tutorial
 *
 * \code{.cpp}
 * real
 * LevelSetFunction( const Matrix< DDRMat > & aPoint )
 * {
 *     return norm( aPoint ) - 1.2;
 * }
 * \endcode
 */
real
LevelSetFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

//------------------------------------------------------------------------------


/*!
 * <h1>Tutorial 2 : Load 2D Mesh from File and refine again"</h1>
 * In the following example, a mesh and a nodal field are loaded
 * from existing files. The refinement is performed and
 * the new mesh is stored into an output.
 *
 * Note that the files created by Tutorial 1 must exist on the file system.
 */
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

//------------------------------------------------------------------------------
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
    Mesh * tInputMesh = tHMR.create_mesh();

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
     * load_matrix_from_binary_file( tInputField->get_node_values(), "Tutorial_1.bin" );
     * \endcode
     */

    // fixme: Make input mesh calculate T-Matrices after initialization.
    //        Then, coefficents can be loaded from binary as well
    load_matrix_from_binary_file( tInputField->get_node_values(), "Tutorial_1.bin" );


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
    Mesh * tOutputMesh = tHMR.create_mesh();

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
     * save_matrix_to_binary_file( tInputField->get_node_values(), "Tutorial_2.bin" );
     * \endcode
     */
    save_matrix_to_binary_file( tInputField->get_node_values(), "Tutorial_2.bin" );

    /*!
     * STK_Implementation is currently not linked to the new MTK interface.
     * This is a workaround for now.
     * \code{.cpp}
     * tHMR.save_to_exodus("Tutorial_2.exo");
     * \endcode
     */
    tHMR.save_to_exodus( "Tutorial_2.exo" );

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

//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}

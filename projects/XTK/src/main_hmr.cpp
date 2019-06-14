/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */
//------------------------------------------------------------------------------

// moris core includes
#include "cl_GE_Core.hpp"
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_HMR_Field.hpp"

//------------------------------------------------------------------------------

// geometry engine
//#include <GEN/src/cl_GEN_Geometry_Engine.hpp>
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

//------------------------------------------------------------------------------
// XTK
#include "xtk/cl_XTK_Model.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "xtk_typedefs.hpp"
#include "cl_Geom_Field.hpp"
#include "cl_Discrete_Level_Set.hpp"

//------------------------------------------------------------------------------



// select namespaces
using namespace moris;
using namespace hmr;
using namespace xtk;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;
/*!
 * \section Tutorial_4: Using the MTK Interface
 *
 * In this example, we will refine a mesh according to the surface of
 * a sphere. When we are done, we will create an MTK interface object.
 *
 */
real
CircleFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}


int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    //------------------------------------------------------------------------------

//    ParameterList tParameters = create_hmr_parameter_list();
//    tParameters.set( "number_of_elements_per_dimension", "2, 2, 2" );


    /*!
     * the domain will have a width of 3 and a height of 3 units
     *
     * \code{.cpp}
     * tParameters.set( "domain_dimensions", "3, 3" );
     * \endcode
     */
//    tParameters.set( "domain_dimensions", "3, 3, 3" );

    /*!
     * set the coordinates of the lower left node on the mesh
     *
     * \code{.cpp}
     * tParameters.set( "domain_offset", "-1.5, -1.5" );
     * \endcode
     */
//    tParameters.set( "domain_offset", "-1.5, -1.5,-1.5" );


    /*!
     * By default, HMR is silent.
     * If this output is set, HMR will print lots of debug information
     * on screen.
     *
     * \code{.cpp}
     * tParameters.set( "verbose", 1 );
     * \endcode
     */
//    tParameters.set( "verbose", 1 );
//    tParameters.set( "truncate_bsplines", 1 );
//
//    tParameters.set( "bspline_orders", "1" );
//    tParameters.set( "lagrange_orders", "1" );
//    tParameters.set( "use_multigrid", 1);

    // order for this example
    //moris::uint tOrder = 1;

    // create parameter object
    moris::hmr::Parameters tParameters;
    tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );
    tParameters.set_domain_dimensions( { { 3 }, { 3 } } );
    tParameters.set_domain_offset( { { 0 }, { 0 } } );
    tParameters.set_verbose( false );
    tParameters.set_multigrid( true );
    tParameters.set_bspline_truncation( true );
    tParameters.set_mesh_orders_simple( 1 );
    tParameters.set_refinement_buffer( 1 );



    //------------------------------------------------------------------------------

    /*!
     * <b> Step 2: HMR object </b>
     */

    /*!
     * All operations such as refining a mesh according to a field,
     * mapping a field onto a new mesh and providing the API of the new mesh
     * are handled by the HMR object.
     *
     * \code{.cpp}
     * HMR tHMR( tParameters );
     * \endcode
     */
    HMR tHMR( tParameters );

    //------------------------------------------------------------------------------

    /*!
     * <b> Step 3: Creating a nodal field and refining according to it</b>
     */

    /*!
     * The following command creates a shared pointer to a field that is
     * called "LevelSet". The datatype is std::shared_ptr<moris::hmr::Field>
     *
     * \code{.cpp}
     * auto tField = tHMR.create_field( "Circle" );
     * \endcode
     */
    auto tField = tHMR.create_field( "Circle" );

    /*!
     * This example uses an analytic level set, which is defined as follows
     *
     * \code{.cpp}
     * real
     * CircleFunction( const Matrix< DDRMat > & aPoint )
     * {
     *     return norm( aPoint ) - 1.2;
     * }
     * \endcode
     *
     * The pointer of this function is passed to the field.
     *
     * \code{.cpp}
     * tField->evaluate_scalar_function( CircleFunction );
     * \endcode
     */
    tField->evaluate_scalar_function( CircleFunction );

    /*!
     * In the next step, we use this field to identify elements that
     * are fully inside the level set, or intersected.
     *
     * This flagging can be repeated with an arbitrary number
     * of fields.
     *
     * \code{.cpp}
     * tHMR.flag_volume_and_surface_elements( tField );
     * \endcode
     */
    tHMR.flag_volume_and_surface_elements( tField );

    /*!
     * One all elements are flagged for the refinement, a procedure is
     * called which performs one refinement step and maps all fields
     * to the new mesh.
     *
     * \code{.cpp}
     * tHMR.perform_refinement_and_map_fields();
     * \endcode
     */
    tHMR.perform_refinement_and_map_fields();

    //------------------------------------------------------------------------------

    /*!
     * <b> Step4: Saving the data</b>
     */

    /*!
     * This command creates an exodus file of the refined mesh.
     *
     * \code{.cpp}
     * tHMR.save_to_exodus( "Mesh1.exo" );
     * \endcode
     */
    tHMR.save_to_exodus( "Mesh1_1.exo" );

    tHMR.save_bsplines_to_vtk("BSplines.vtk");

    /*!
     * One can also recall the state of the mesh before the last
     * refinement, and save it into a file.
     *
     * \code{.cpp}
     * tHMR.save_last_step_to_exodus( "LastStep.exo" );
     * \endcode
     */
    tHMR.save_last_step_to_exodus( "LastStep.exo" );

    /*!
     * The database that defines the mesh in its refined state saved using this command.
     * Note that field information is not stored.
     *
     * \code{.cpp}
     * tHMR.save_to_hdf5( "Database.hdf5" );
     * \endcode
     */
    tHMR.save_to_hdf5( "Database.hdf5" );

    /*!
     * For third party applications, the interpolation coefficients
     * can be saved into an output file as well
     * \code{.cpp}
     * tHMR.save_coeffs_to_hdf5_file( "TMatrix.hdf5" );
     * \endcode
     */
    tHMR.save_coeffs_to_hdf5_file( "TMatrix.hdf5" );

    /*!
     * Each individual field can be stored into an hdf5 file on its own.
     * These files can be exchanged with other codes such as MATLAB
     * \code{.cpp}
     * tField->save_field_to_hdf5("Circle.hdf5");
     * \endcode
     */
    tField->save_field_to_hdf5( "Circle.hdf5" );
    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( );

//    for(moris::uint  i = 0; i < tMesh->get_num_entities(EntityRank::NODE); i++)
//    {
//        std::cout<<" i = "<<i<<std::endl;
//
//        moris::print(tMesh->get_bspline_inds_of_node_loc_ind(i,EntityRank::BSPLINE_1),"bspline inds");
//        moris::print(tMesh->get_t_matrix_of_node_loc_ind(i,EntityRank::BSPLINE_1),"bspline t matrix");
//    }
//
//
//    moris::moris_id tElementInd = 8;
//    Matrix<IndexMat> tFaces = tMesh->get_entity_connected_to_entity_loc_inds(tElementInd,EntityRank::ELEMENT,EntityRank::FACE);
//    Matrix<IndexMat> tNeighborElements = tMesh->get_elements_connected_to_element_and_face_ind_loc_inds(tElementInd);
//    moris::print(tNeighborElements,"tNeighborElements");
//    moris::print(tFaces,"my faces");
//    for(moris::uint i =0 ; i <tNeighborElements.n_cols(); i++)
//    {
//        moris::print(tMesh->get_entity_connected_to_entity_loc_inds(tNeighborElements(0,i),EntityRank::ELEMENT,EntityRank::FACE),"neighbors faces");
//    }

    xtk::Geom_Field tFieldAsGeom(tField);

//    moris::mtk::Mesh* tMeshData   = moris::mtk::create_mesh( MeshType::STK, "Mesh1.exo" );
//    std::string tLSFName            = "Circle";
//    xtk::Discrete_Level_Set tLevelSetMesh(tMeshData,{tLSFName});

    // Tell the geometry engine about the discrete field mesh and how to interpret phases
    Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine tGeometryEngine(tFieldAsGeom,tPhaseTable);


    // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tMesh.get(),tGeometryEngine);
    tXTKModel.mSameMesh = true;
    tXTKModel.mVerbose = true;

    // Do the cutting
    tXTKModel.decompose(tDecompositionMethods);

    tXTKModel.unzip_interface();

//    tXTKModel.perform_basis_enrichment();


    Output_Options tOutputOptions;
    tOutputOptions.mAddNodeSets = false;
    tOutputOptions.mAddSideSets = false;

    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);


    std::string tPrefix = std::getenv("MORISOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/xtk_hmr_output.e";
    tCutMeshData->create_output_mesh(tMeshOutputFile);
    delete tCutMeshData;
//    delete tMeshData;

    //------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}


/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */
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
#include "fn_print.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

//------------------------------------------------------------------------------
// XTK
#include "cl_XTK_Model.hpp"
#include "xtk_typedefs.hpp"

//------------------------------------------------------------------------------
#include "cl_GEN_User_Defined_Geometry.hpp"
#include "cl_GEN_Phase_Table.hpp"

#include "fn_PRM_HMR_Parameters.hpp"


// select namespaces
using namespace moris;
using namespace hmr;
using namespace xtk;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;


moris::real
LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{

    real mXn = 0;
    real mYn = 0;
    real mZn = 1.0;
    real mXc = 1.0;
    real mYc = 1.0;
    real mZc = 3.7;
    return mXn*(aPoint(0)-mXc) + mYn*(aPoint(1)-mYc) + mZn*(aPoint(2)-mZc);
}

moris::real LevelSetFunctionGeometry(
        const moris::Matrix<moris::DDRMat>& aCoordinates,
        const moris::Cell<moris::real*>& aParameters)
{
    return LevelSetFunction(aCoordinates);
}

moris::real
LevelSetFunction2( const moris::Matrix< moris::DDRMat > & aPoint )
{
    return norm( aPoint ) - 2.611;
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

    //    moris::uint tBplineOrder = 1;
    moris::uint tLagrangeOrder = 1;
    //    moris::uint tMyCoeff = 1;

    moris::ParameterList tParameters = prm::create_hmr_parameter_list();

    tParameters.set( "number_of_elements_per_dimension", "1, 1, 4" );
    tParameters.set( "domain_dimensions", "1, 1, 4" );
    tParameters.set( "domain_offset", "0.0, 0.0, 0.0" );
    tParameters.set( "domain_sidesets", "5");
    tParameters.set( "verbose", 0 );
    tParameters.set( "truncate_bsplines", 1 );
    tParameters.set( "bspline_orders", "1" );
    tParameters.set( "lagrange_orders", "1" );
    tParameters.set( "use_multigrid", 0 );
    tParameters.set( "refinement_buffer", 2 );
    tParameters.set( "staircase_buffer", 2 );

    hmr::HMR tHMR( tParameters );

    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

    // create field
    std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle1", tLagrangeOrder );
    //    std::shared_ptr< moris::hmr::Field > tField2 = tMesh->create_field( "Circle2", tLagrangeOrder );

    tField->evaluate_scalar_function( LevelSetFunction );
    //    tField2->evaluate_scalar_function( LevelSetFunction2 );

    for( uint k=0; k<3; ++k )
    {
        tHMR.flag_surface_elements_on_working_pattern( tField );
        tHMR.perform_refinement_based_on_working_pattern( 0 );

        tField->evaluate_scalar_function( LevelSetFunction );
    }

    tHMR.finalize();

    Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh(
            tLagrangeOrder, tHMR.mParameters->get_lagrange_output_pattern() );

    std::cout<<"Num Nodes ="<<tMesh->get_num_nodes()<<std::endl;
    std::cout<<"Num Cells ="<<tMesh->get_num_elems()<<std::endl;

    moris::Cell<std::shared_ptr<moris::ge::Geometry>> tGeometryVector(1);
    tGeometryVector(0) = std::make_shared<ge::User_Defined_Geometry>(Matrix<DDRMat>(0, 0), &(LevelSetFunctionGeometry));

    // Tell the geometry engine about the discrete field mesh and how to interpret phases
    moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
    tGeometryEngineParameters.mGeometries = tGeometryVector;
    moris::ge::Geometry_Engine tGeometryEngine(tInterpMesh, tGeometryEngineParameters);

    // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tInterpMesh,&tGeometryEngine);
    tXTKModel.mVerbose  =  false;

    // Do the cutting
    tXTKModel.decompose(tDecompositionMethods);

    //    print_cut_mesh(tXTKModel);

    //    tXTKModel.unzip_interface();
    //
    //    tXTKModel.perform_basis_enrichment();

    Output_Options tOutputOptions;
    tOutputOptions.mAddNodeSets = false;
    tOutputOptions.mAddSideSets = true;
    tOutputOptions.mAddClusters = true;

    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

    std::string tMeshOutputFile = "./xtk_exo/xtk_hmr_output.e";
    tCutMeshData->create_output_mesh(tMeshOutputFile);
    delete tCutMeshData;
    delete tInterpMesh;
    //    delete tMeshData;

    //------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}


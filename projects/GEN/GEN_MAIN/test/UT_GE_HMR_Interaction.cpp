#include <catch.hpp>

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
//#include "cl_Geom_Field.hpp"
#include "typedefs.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src


#include "cl_MTK_Mesh_Manager.hpp"

#include "HDF5_Tools.hpp"

// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Intersection_Object_Line.hpp"
#include "cl_GE_Node.hpp"

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Geom_Field.hpp"

// LINALG includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Vertex.hpp"

#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

using namespace moris;
using namespace hmr;

namespace ge
{

moris::real
LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint, const moris::Cell< moris::real > aConst )
{
    return norm( aPoint ) - 0.9;
}

TEST_CASE("GE_HMR_Interaction_00","[moris],[GE],[GE_HMR_Interaction]")
{
    if(par_size() == 1)
    {
        for( moris::uint tOrder=1; tOrder<=1; tOrder++ )
        {
            uint tLagrangeMeshIndex = 0;
            // empty container for B-Spline meshes
            moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { {4}, {4} } );

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders  ( { {tOrder} });
            tParameters.set_lagrange_patterns({ {2} });

            tParameters.set_bspline_orders   ( { {tOrder}, {tOrder}} );
            tParameters.set_bspline_patterns ( { {0}, {1}} );

            tParameters.set_staircase_buffer( tOrder );

            tParameters.set_initial_refinement( 2 );

            Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { {0}, {1} };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // manually select output pattern
            tDatabase->set_activation_pattern( 0 );

            tHMR.perform_initial_refinement( 0 );

//            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_0_initial.vtk");

            // manually select output pattern
            tDatabase->set_activation_pattern( 1 );

            // refine the last element three times
            // fixme: change this to 2
            for( uint tLevel = 0; tLevel < 1; ++tLevel )
            {
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }

//            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_1_initial.vtk");

            tDatabase->unite_patterns( 0, 1, 2 );

//            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_2_initial.vtk");

            tDatabase->update_bspline_meshes();

            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
            std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 2,*tInterpMesh );

            mtk::Mesh_Manager tMesh;
            uint tMeshIndex = tMesh.register_mesh_pair( tInterpMesh.get(), tIntegrationMesh.get() );

            //--------------------------------------------------------------------------------------------------

            // input parameters for the circle LS
//            moris::Cell< real > tCircleInputs = { {0.9}, {0}, {0} };
            //------------------------------------------------------------------------------

            moris::ge::Ge_Factory tFactory;
            std::shared_ptr< moris::ge::Geometry > tGeom = tFactory.set_geometry_type( moris::ge::GeomType::ANALYTIC );

            tGeom->set_my_mesh( &tMesh );

//            moris_index tSubIndex = tGeom->set_analytical_function( LevelSetFunction, tCircleInputs );
            moris_index tSubIndex = tGeom->set_analytical_function( LevelSetFunction );

            moris::ge::GE_Core tGeometryEngine;
            moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom );

            uint tNumOfIPNodes = tIntegrationMesh->get_num_nodes();

            Matrix< DDRMat > tFieldData = tGeometryEngine.get_field_vals( tMyGeomIndex, tSubIndex );
//            print(tFieldData, "tFieldData");

            tHMR.based_on_field_put_elements_on_queue( tFieldData, tLagrangeMeshIndex );

            tDatabase->get_background_mesh()->perform_refinement( 1 );


//            tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE, 1 );
            tHMR.update_refinement_pattern( 1 );

            moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes1;

            // create factory
            moris::hmr::Factory tFactory_HMR;

            moris::hmr::Lagrange_Mesh_Base* tLagrangeMesh =  tFactory_HMR.create_lagrange_mesh( &tParameters,
                                                                                                tDatabase->get_background_mesh(),
                                                                                                 tBSplineMeshes1,
                                                                                                 1,
                                                                                                 1 );

             // output to exodus
//             STK * tSTK = tLagrangeMesh->create_stk_object(0);
//             tSTK->save_to_file( "GE_HMR_Mesh.g");
//             delete tSTK;

             REQUIRE( tLagrangeMesh->get_number_of_nodes_on_proc()  == 59 );

        delete tLagrangeMesh;
    }
    }
}
//------------------------------------------------------------------------------

TEST_CASE( "GE_HMR_Interaction_01","[GE_HMR_Interaction_Gyroid]" )
{
/*
    if(par_size() == 1)
    {
        uint tLagrangeMeshIndex = 0;
        // empty container for B-Spline meshes
        moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

        // create settings object
        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {10}, {10}, {10} } );
        tParameters.set_domain_dimensions( 10,10,10 );
        tParameters.set_domain_offset( -5, -5, -5 );

        tParameters.set_bspline_truncation( true );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns( { {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_output_meshes( { {0} } );

        tParameters.set_staircase_buffer( 1 );

        tParameters.set_initial_refinement( 1 );

        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        // create the HMR object by passing the settings to the constructor
        moris::hmr::HMR tHMR( tParameters );

        // std::shared_ptr< Database >
        auto tDatabase = tHMR.get_database();

        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );   // HMR Lagrange mesh
        //==============================
        std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "gyroid", tLagrangeMeshIndex);

        tField->evaluate_scalar_function( moris::ge::getDistanceToGyroidsMassive );

        for( uint k=0; k<2; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( moris::ge::getDistanceToGyroidsMassive );
        }
            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_initial_3x3x3.vtk");
        //==============================
        tDatabase->update_bspline_meshes();
        tDatabase->update_lagrange_meshes();

        // calculate T-Matrices etc
        tField->evaluate_scalar_function( moris::ge::getDistanceToGyroidsMassive );
        tDatabase->finalize();
//==============================
        tHMR.save_to_exodus( 0, "gyroid_general_geomEng.g" );

        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        moris::ge::GEN_Geom_Field tFieldAsGeom(tField);

        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};

        size_t tModelDimension = 3;
        moris::ge::GEN_Phase_Table  tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::GEN_Geometry_Engine  tGeometryEngine( tGeometryVector,tPhaseTable,tModelDimension );
        xtk::Model                  tXTKModel( tModelDimension,tInterpMesh.get(),tGeometryEngine );
        tXTKModel.mVerbose = false;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        //tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

        //++++++++++++++++++++++++++++++++
//        xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

        // Declare the fields related to enrichment strategy in output options
//        Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();
        //++++++++++++++++++++++++++++++++

               // output solution and meshes
                       xtk::Output_Options tOutputOptions1;
                       tOutputOptions1.mAddNodeSets = false;
                       tOutputOptions1.mAddSideSets = false;
                       tOutputOptions1.mAddClusters = false;

//                       tOutputOptions1.mRealElementExternalFieldNames = tEnrichmentFieldNames;

                       //++++++++++++++++++++++++++++++++
                       // add solution field to integration mesh
                       std::string tIntegSolFieldName = "solution";
                       tOutputOptions1.mRealNodeExternalFieldNames = {tIntegSolFieldName};
                       //++++++++++++++++++++++++++++++++

                       moris::mtk::Integration_Mesh* tIntegMesh11 = tXTKModel.get_output_mesh(tOutputOptions1);

                       //++++++++++++++++++++++++++++++++
//                       tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tIntegMesh11);
                       //++++++++++++++++++++++++++++++++

//                       for(moris::uint i = 0; i < tIntegMesh11->get_num_entities(EntityRank::NODE); i++)
//                       {
//                           moris::moris_id tID = tIntegMesh11->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//                       }


//                       std::string tMeshOutputFile1 = "./output_general_geomEng.e";
//                       tIntegMesh11->create_output_mesh(tMeshOutputFile1);

//       xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//       xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//==============================
       delete tIntegMesh11;

    }
*/
}

}   // ge namepsace

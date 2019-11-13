/*
 * UT_XTK_HMR_2D.cpp
 *
 *  Created on: Sep 10, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"

//#include "cl_Geom_Field.hpp"
//#include "cl_Plane.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src



#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "../projects/GEN/src/geometry/cl_GEN_Geom_Field.hpp"
#include "../projects/GEN/src/geometry/cl_GEN_Geometry.hpp"
#include "../projects/GEN/src/geometry/cl_GEN_Plane.hpp"

#include "fn_norm.hpp"

namespace xtk
{
moris::real
CircleFuncXTKHMR2D(const moris::Matrix< moris::DDRMat > & aPoint )
{

    moris::real mXCenter = 0;
    moris::real mYCenter = 0;
    moris::real mRadius = 1.1;

    return    (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
            + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
            - (mRadius * mRadius);
}

TEST_CASE("2D XTK WITH HMR","[XTK_HMR_2D]")
{
    if(par_size()<=1)
    {
        std::string tFieldName = "Cylinder";

         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         moris::hmr::Parameters tParameters;

         tParameters.set_number_of_elements_per_dimension( { {2}, {2}} );
         tParameters.set_domain_dimensions({ {2}, {2} });
         tParameters.set_domain_offset({ {-1.0}, {-1.0} });
         tParameters.set_bspline_truncation( true );

         tParameters.set_output_meshes( { {0} } );

         tParameters.set_lagrange_orders  ( { {1} });
         tParameters.set_lagrange_patterns({ {0} });

         tParameters.set_bspline_orders   ( { {1} } );
         tParameters.set_bspline_patterns ( { {0} } );

         tParameters.set_side_sets({{1},{2},{3},{4} });

         tParameters.set_union_pattern( 2 );
         tParameters.set_working_pattern( 3 );

         tParameters.set_refinement_buffer( 2 );
         tParameters.set_staircase_buffer( 2);

         Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
         tLagrangeToBSplineMesh( 0 ) = { {0} };

         tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

         hmr::HMR tHMR( tParameters );

         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

         // create field
         std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

         tField->evaluate_scalar_function( CircleFuncXTKHMR2D );

         for( uint k=0; k<3; ++k )
         {
             tHMR.flag_surface_elements_on_working_pattern( tField );
             tHMR.perform_refinement_based_on_working_pattern( 0 );

             tField->evaluate_scalar_function( CircleFuncXTKHMR2D );
         }

         tHMR.finalize();

         tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_ip.e" );

         std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

         moris::ge::GEN_Geom_Field tFieldAsGeom(tField);

         moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};

         size_t tModelDimension = 2;
         moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
         moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
         Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
         tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

//        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

        // output to exodus file ----------------------------------------------------------
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        std::string tMeshOutputFile ="./xtk_exo/xtk_hmr_2d_ig.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);

        delete tCutMeshData;
    }
}

TEST_CASE("2D XTK WITH HMR WEIRD INTERSECTION","[XTK_HMR_2D_WI]")
{
    if(par_size()<=1)
    {
        std::string tFieldName = "Cylinder";

         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         moris::hmr::Parameters tParameters;

         tParameters.set_number_of_elements_per_dimension( { {2}, {2}} );
         tParameters.set_domain_dimensions({ {2}, {2} });
         tParameters.set_domain_offset({ {-1.0}, {-1.0} });
         tParameters.set_bspline_truncation( true );

         tParameters.set_output_meshes( { {0} } );

         tParameters.set_lagrange_orders  ( { {1} });
         tParameters.set_lagrange_patterns({ {0} });

         tParameters.set_bspline_orders   ( { {1} } );
         tParameters.set_bspline_patterns ( { {0} } );

         tParameters.set_side_sets({{1},{2},{3},{4} });

         tParameters.set_union_pattern( 2 );
         tParameters.set_working_pattern( 3 );

         tParameters.set_refinement_buffer( 2 );
         tParameters.set_staircase_buffer( 2);

         Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
         tLagrangeToBSplineMesh( 0 ) = { {0} };

         tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

         hmr::HMR tHMR( tParameters );

         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

         // create field
         std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

         tField->evaluate_scalar_function( CircleFuncXTKHMR2D );

         for( uint k=0; k<2; ++k )
         {
             tHMR.flag_surface_elements_on_working_pattern( tField );
             tHMR.perform_refinement_based_on_working_pattern( 0 );

             tField->evaluate_scalar_function( CircleFuncXTKHMR2D );
         }

         tHMR.finalize();

         tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_wi_2d_ip.e" );

         std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

         // create a plane which intentionally intersects from fine to coarse
         moris::Matrix<moris::DDRMat> tCenters = {{ 0.1,0.1 }};
         moris::Matrix<moris::DDRMat> tNormals = {{ 1.0,0.0 }};
         moris::ge::Plane<2> tPlane(tCenters,tNormals);

         moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tPlane};

         size_t tModelDimension = 2;
         moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
         moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
         Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
         tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

        // output to exodus file ----------------------------------------------------------
        xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

         // Declare the fields related to enrichment strategy in output options
         Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();

        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tIntegMesh1);

        std::string tMeshOutputFile ="./xtk_exo/xtk_hmr_2d_wi_ig.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tIntegMesh1;
    }
}

}

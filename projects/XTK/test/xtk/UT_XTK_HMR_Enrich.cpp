/*
 * UT_XTK_HMR_Enrich.cpp
 *
 *  Created on: Sep 30, 2019
 *      Author: doble
 */




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
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

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

#include "cl_GEN_Geom_Field.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Plane.hpp"

#include "fn_norm.hpp"

namespace xtk
{
moris::real
MultiCircle(const moris::Matrix< moris::DDRMat > & aPoint )
{

    moris::real mXCenter = 0.3333;
    moris::real mYCenter = 0.3333;
    moris::real mRadius = 0.22;


    real val1 =   (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
                    + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
                    - (mRadius * mRadius);

    mXCenter = -0.3333;
    mYCenter = -0.3333;

    real val2 = (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
                            + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
                            - (mRadius * mRadius);


    moris::real mXC =  0.422;
    moris::real mYC = -0.422;
    moris::real mNx = 1.0;
    moris::real mNy = -1.0;
    moris::real val3 = -(mNx*(aPoint(0)-mXC) + mNy*(aPoint(1)-mYC));

    mXC = -0.422;
    mYC = 0.422;
    mNx = 1.0;
    mNy = -1.0;
    moris::real val4 = (mNx*(aPoint(0)-mXC) + mNy*(aPoint(1)-mYC));

    return std::min(val1,std::min(val2,std::min(val3,val4)));

}

TEST_CASE("2D XTK WITH HMR No truncation enrichment","[XTK_HMR_ENR_2D]")
{
    if(par_size()<=1)
    {
        std::string tFieldName = "Cylinder";

         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         moris::hmr::Parameters tParameters;

         tParameters.set_number_of_elements_per_dimension( { {3}, {3}} );
         tParameters.set_domain_dimensions({ {2}, {2} });
         tParameters.set_domain_offset({ {-1.0}, {-1.0} });
         tParameters.set_bspline_truncation( true );

         tParameters.set_output_meshes( { {0} } );

         tParameters.set_lagrange_orders  ( { {2} });
         tParameters.set_lagrange_patterns({ {0} });

         tParameters.set_bspline_orders   ( { {3} } );
         tParameters.set_bspline_patterns ( { {0} } );

         tParameters.set_side_sets({{1},{2},{3},{4} });

         tParameters.set_union_pattern( 2 );
         tParameters.set_working_pattern( 3 );

         tParameters.set_refinement_buffer( 1 );
         tParameters.set_staircase_buffer( 1 );

         Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
         tLagrangeToBSplineMesh( 0 ) = { {0} };

         tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

         hmr::HMR tHMR( tParameters );

         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

         // create field
         std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

         tField->evaluate_scalar_function( MultiCircle );

         for( uint k=0; k<2; ++k )
         {
             tHMR.flag_surface_elements_on_working_pattern( tField );
             tHMR.perform_refinement_based_on_working_pattern( 0 );
             tField->evaluate_scalar_function( MultiCircle );
         }

         tHMR.finalize();

         tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_enr_ip.e" );

         std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

         moris::ge::GEN_Geom_Field tFieldAsGeom(tField);

         moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};

         size_t tModelDimension = 2;
         moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
         moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
         Model tXTKModel(tModelDimension,tInterpMesh.get(),&tGeometryEngine);
         tXTKModel.mVerbose  =  false;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);

        Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

        // Declare the fields related to enrichment strategy in output options
        Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tCutMeshData);

        std::string tMeshOutputFile ="./xtk_exo/xtk_hmr_2d_enr_trun_ig.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);

        delete tCutMeshData;
    }
}


moris::real
CircleMultiMat(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXCenter = 0.01;
    moris::real mYCenter = 0.01;
    moris::real mRadius = 0.61;


    return  (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
                    + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
                    - (mRadius * mRadius);
}

moris::real
PlaneMultiMat(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXC = 0.01;
    moris::real mYC = 0.01;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return (mNx*(aPoint(0)-mXC) + mNy*(aPoint(1)-mYC));
}

TEST_CASE("2D XTK WITH HMR Multi-Mat","[XTK_HMR_MULTI_2D]")
{
    if(par_size()<=1)
    {
        std::string tFieldName = "Geometry";

         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         moris::hmr::Parameters tParameters;

         tParameters.set_number_of_elements_per_dimension( { {3}, {1}} );
         tParameters.set_domain_dimensions({ {6}, {2} });
         tParameters.set_domain_offset({ {-3.0}, {-1.0} });
         tParameters.set_bspline_truncation( true );

         tParameters.set_output_meshes( { {0} } );

         tParameters.set_lagrange_orders  ( { {2} });
         tParameters.set_lagrange_patterns({ {0} });

         tParameters.set_bspline_orders   ( { {3} } );
         tParameters.set_bspline_patterns ( { {0} } );

         tParameters.set_side_sets({{1},{2},{3},{4} });
         tParameters.set_max_refinement_level( 2 );
         tParameters.set_union_pattern( 2 );
         tParameters.set_working_pattern( 3 );

         tParameters.set_refinement_buffer( 2 );
         tParameters.set_staircase_buffer( 2 );

         Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
         tLagrangeToBSplineMesh( 0 ) = { {0} };

         tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

         hmr::HMR tHMR( tParameters );

         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

         // create field
         std::shared_ptr< moris::hmr::Field > tPlaneField  = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

         tPlaneField->evaluate_scalar_function( PlaneMultiMat );

         for( uint k=0; k<2; ++k )
         {
             tHMR.flag_surface_elements_on_working_pattern( tPlaneField );
             tHMR.perform_refinement_based_on_working_pattern( 0 );
             tPlaneField->evaluate_scalar_function( PlaneMultiMat );
         }

         std::shared_ptr< moris::hmr::Field > tCircleField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
         tCircleField->evaluate_scalar_function( CircleMultiMat );
         for( uint k=0; k<2; ++k )
         {
             tHMR.flag_surface_elements_on_working_pattern( tCircleField );
             tHMR.perform_refinement_based_on_working_pattern( 0 );
             tCircleField->evaluate_scalar_function( CircleMultiMat );

         }

         tPlaneField->evaluate_scalar_function( PlaneMultiMat );
         tCircleField->evaluate_scalar_function( CircleMultiMat );
         tHMR.finalize();

         tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_enr_ip2.e" );

         std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

         moris::ge::GEN_Geom_Field tPlaneFieldAsGeom(tPlaneField);
         moris::ge::GEN_Geom_Field tCircleFieldAsGeom(tCircleField);

         moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tCircleFieldAsGeom,&tPlaneFieldAsGeom};

         size_t tModelDimension = 2;
         moris::ge::GEN_Phase_Table tPhaseTable (2,  Phase_Table_Structure::EXP_BASE_2);
         moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
         Model tXTKModel(tModelDimension,tInterpMesh.get(),&tGeometryEngine);
         tXTKModel.mVerbose  =  false;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);

        Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

        Enriched_Integration_Mesh & tEnrInteg = tXTKModel.get_enriched_integ_mesh(0);

        // Declare the fields related to enrichment strategy in output options
        Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tCutMeshData);

        std::string tMeshOutputFile ="./xtk_exo/xtk_hmr_2d_enr_trun_ig.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);

        delete tCutMeshData;
    }
}
}


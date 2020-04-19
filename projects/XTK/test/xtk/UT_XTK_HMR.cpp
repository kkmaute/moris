/*
 * UT_XTK_HMR.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
//#include "cl_Geom_Field.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Checker.hpp"

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

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

#include "cl_GEN_Geometry_Discrete.hpp"
#include "cl_GEN_Geometry_Field_HMR.hpp"

#include "fn_norm.hpp"


namespace moris
{


moris::real
LevelSetSphereCylinder(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::Matrix<moris::DDRMat> aCenter = {{0.0},{0.0},{0.0}};
    moris::Matrix<moris::DDRMat> aAxis   = {{0.0},{1.0},{0.0}};
    moris::real aRad = 0.77;
    moris::real aLength = 5;

    MORIS_ASSERT(aCenter.numel() == 3,"Centers need to have length 3");
    MORIS_ASSERT(aAxis.numel() == 3, "axis need to have length 3");

    Cell<moris::real> relativePosition = {(aPoint(0) - aCenter(0)),(aPoint(1) - aCenter(1)),(aPoint(2) - aCenter(2))};
    moris::real lsFromLeft = (relativePosition(0)*(-aAxis(0)) + relativePosition(1)*(-aAxis(1))+ relativePosition(2)*(-aAxis(2))) - aLength/2.0;
    moris::real lsFromRight = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2))) - aLength/2.0;

    moris::real axialCrd = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2)));
    Cell<moris::real> radDir = {(relativePosition(0) - aAxis(0)*axialCrd), (relativePosition(1) - aAxis(1)*axialCrd),(relativePosition(2) - aAxis(2)*axialCrd)};
    moris::real radDist = std::pow(radDir(0)*radDir(0)+radDir(1)*radDir(1)+radDir(2)*radDir(2), 0.5);
    moris::real lsFromRad = radDist - aRad;

    return -std::max(std::max(lsFromLeft, lsFromRight), lsFromRad);
}

moris::real
LevelSetPlaneFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{

    real mXn = 0.0;
    real mYn = 0.0;
    real mZn = 1.0;
    real mXc = 0.0;
    real mYc = 0.0;
    real mZc = 1.8;
    return mXn*(aPoint(0)-mXc) + mYn*(aPoint(1)-mYc) + mZn*(aPoint(2)-mZc);
}

TEST_CASE("XTK HMR Test","[XTK_HMR]")
{
//    if(par_size() == 1)
//    {
//        std::string tFieldName = "Cylinder";
//
//        moris::uint tLagrangeMeshIndex = 0;
//        moris::uint tBSplineMeshIndex = 0;
//
//        moris::hmr::Parameters tParameters;
//
//        tParameters.set_number_of_elements_per_dimension( { {4}, {4}, {4} } );
//        tParameters.set_domain_dimensions({ {2}, {2}, {4} });
//        tParameters.set_domain_offset({ {-1.0}, {-1.0}, {-2.0} });
//        tParameters.set_side_sets({ {5}, {6} });
//
//
//        tParameters.set_bspline_truncation( true );
//
//        tParameters.set_lagrange_orders  ( { {1} });
//        tParameters.set_lagrange_patterns( { {0} });
//
//        tParameters.set_bspline_orders   ( { {1} } );
//        tParameters.set_bspline_patterns ( { {0} } );
//
//        tParameters.set_output_meshes( { {0} } );
//        //        tParameters.set_lagrange_input_mesh( { { 0 } } );
//
//        tParameters.set_staircase_buffer( 1 );
//
//        tParameters.set_initial_refinement( 0 );
//
//        tParameters.set_number_aura( true );
//
//        Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
//        tLagrangeToBSplineMesh( 0 ) = { {0} };
//
//        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );
//
//        // create the HMR object by passing the settings to the constructor
//        moris::hmr::HMR tHMR( tParameters );
//
//        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//        // create field
//        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
//
//        tField->evaluate_scalar_function( LevelSetSphereCylinder );
//
//        for( uint k=0; k<3; ++k )
//        {
//            tHMR.flag_surface_elements_on_working_pattern( tField );
//            tHMR.perform_refinement_based_on_working_pattern(0 );
//
//            tField->evaluate_scalar_function( LevelSetSphereCylinder );
//        }
//
//        tHMR.finalize();
//
//        tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_ghost_interp.e" );
//
//        hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//        hmr::Lagrange_Mesh_Base * tLMB = tInterpMesh->get_lagrange_mesh();
//
//        for(moris::uint  i = 0 ; i < tLMB->mFacets.size(); i++)
//        {
//            if(! (tLMB->mFacets(i) == nullptr) )
//            {
//                hmr::Element * tMaster = tLMB->mFacets(i)->get_hmr_master();
//                hmr::Element * tSlave  = tLMB->mFacets(i)->get_hmr_slave();
//
//                if(! (tMaster == nullptr) && ! (tSlave == nullptr)  )
//                {
//                    moris_id tMasterId = tMaster->get_index();
//                    moris_id tSlaveId  = tSlave->get_index();
//
//                    if(tMasterId != MORIS_INDEX_MAX)
//                    {
//                        tMasterId = tMasterId + 1;
//                    }
//
//                    if(tSlaveId != MORIS_INDEX_MAX)
//                    {
//                        tSlaveId = tSlaveId + 1;
//                    }
//                }
//            }
//        }
//
//
//
//
//
//
//
//        moris::ge::Geometry_Field_HMR tFieldAsGeom(tField);
//
//        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//        // Tell the geometry engine about the discrete field mesh and how to interpret phases
//        moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//        moris::ge::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);
//
//        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
//        size_t tModelDimension = 3;
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
//        xtk::Model tXTKModel(tModelDimension,tInterpMesh,&tGeometryEngine);
//        tXTKModel.mSameMesh = true;
//        tXTKModel.mVerbose  =  false;
//
//        // Do the cutting
//        tXTKModel.decompose(tDecompositionMethods);
//
//        // Perform the enrichment
//        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
//
//        // perform ghost stabilization
////        tXTKModel.construct_face_oriented_ghost_penalization_cells();
//
//
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = true;
//        tOutputOptions.mAddClusters = false;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldName = "solution";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//
//        std::string tOutputFile = "./xtk_exo/xtk_hmr_cut.exo";
//        tIntegMesh1->create_output_mesh(tOutputFile);
//
//        delete tIntegMesh1;
//        delete tInterpMesh;
//    }
}
}


/*
 * UT_XTK_HMR.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Vertex_Enrichment.hpp"
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

#include "cl_GEN_User_Defined_Geometry.hpp"

#include "fn_norm.hpp"

#include <map>
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

real LevelSetSphereCylinderGeometry(const Matrix<DDRMat>& aCoordinates, const Cell<real*>& aParameters)
{
    return LevelSetSphereCylinder(aCoordinates);
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

        //iterate through linear quadratic cubic
        for( moris::uint iOrder = 1; iOrder < 2; iOrder ++)
        {

            std::string tFieldName = "Cylinder";

            moris::uint tLagrangeMeshIndex = 0;
            moris::uint tBSplineMeshIndex = 0;

            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { {6}, {6}, {6} } );
            tParameters.set_domain_dimensions({ {2}, {2}, {4} });
            tParameters.set_domain_offset({ {-1.0}, {-1.0}, {-2.0} });
            tParameters.set_side_sets({ {5}, {6} });


            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders  ( { {iOrder} });
            tParameters.set_lagrange_patterns( { {0} });

            tParameters.set_bspline_orders   ( { {iOrder} } );
            tParameters.set_bspline_patterns ( { {0} } );

            tParameters.set_output_meshes( { {0} } );
            //        tParameters.set_lagrange_input_mesh( { { 0 } } );

            tParameters.set_staircase_buffer( 1 );

            tParameters.set_initial_refinement( 0 );

            tParameters.set_number_aura( true );

            Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { {0} };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LevelSetSphereCylinder );

            for( uint k=0; k<1; ++k )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern(0 );

                tField->evaluate_scalar_function( LevelSetSphereCylinder );
            }

            tHMR.finalize();

            hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

            hmr::Lagrange_Mesh_Base * tLMB = tInterpMesh->get_lagrange_mesh();

            for(moris::uint  i = 0 ; i < tLMB->mFacets.size(); i++)
            {
                if(! (tLMB->mFacets(i) == nullptr) )
                {
                    hmr::Element * tMaster = tLMB->mFacets(i)->get_hmr_master();
                    hmr::Element * tSlave  = tLMB->mFacets(i)->get_hmr_slave();

                    if(! (tMaster == nullptr) && ! (tSlave == nullptr)  )
                    {
                        moris_id tMasterId = tMaster->get_index();
                        moris_id tSlaveId  = tSlave->get_index();

                        if(tMasterId != MORIS_INDEX_MAX)
                        {
                            tMasterId = tMasterId + 1;
                        }

                        if(tSlaveId != MORIS_INDEX_MAX)
                        {
                            tSlaveId = tSlaveId + 1;
                        }
                    }
                }
            }

            Cell< std::shared_ptr<ge::Geometry> > tGeometryVector(1);
            tGeometryVector(0) = std::make_shared<moris::ge::User_Defined_Geometry>(Matrix<DDRMat>(0, 0), &(LevelSetSphereCylinderGeometry));

            size_t tModelDimension = 3;
            moris::ge::Phase_Table tPhaseTable (1, moris::ge::Phase_Table_Structure::EXP_BASE_2);
            moris::ge::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
            xtk::Model tXTKModel(tModelDimension,tInterpMesh,&tGeometryEngine);
            tXTKModel.mVerbose  =  true;

            // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};

            // Do the cutting
            tXTKModel.decompose(tDecompositionMethods);

            // Perform the enrichment
            tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);

            xtk::Enriched_Interpolation_Mesh & tEnrIpMesh  = tXTKModel.get_enriched_interp_mesh(0);

            Matrix<IndexMat> const & tEnrCoeffs = tEnrIpMesh.get_enriched_coefficient_local_to_global_map(0);

            // perform ghost stabilization
            tXTKModel.construct_face_oriented_ghost_penalization_cells();


            // get the enriched basis coeffs

            // Check for double ids/indices
            std::unordered_map<moris_index,moris_index> tBasisIdMap;
            for(moris::uint i = 0 ; i < tEnrCoeffs.numel(); i++)
            {
                MORIS_ASSERT(tBasisIdMap.find(tEnrCoeffs(i)) == tBasisIdMap.end(),"Id already in the map");

                tBasisIdMap[tEnrCoeffs(i)] = i;

            }

            moris::uint tNumVerts = tEnrIpMesh.get_num_entities(EntityRank::NODE,0);


            Cell<moris_index> tOwnedIndices;
            Cell<moris_index> tOwnedIds;

            // iterate through vertices
            for(moris::uint iV = 0; iV < tNumVerts; iV++)
            {
                // get the vertex
                xtk::Interpolation_Vertex_Unzipped & tVertex = tEnrIpMesh.get_xtk_interp_vertex((moris_index)iV);
                if(tVertex.get_base_vertex()->has_interpolation(1))
                {
                    xtk::Vertex_Enrichment* tInterp =  tVertex.get_xtk_interpolation(0);

                    Matrix<IndexMat> tOwners  = tInterp->get_owners();
                    Matrix<IndexMat> tIndices = tInterp->get_indices();
                    Matrix<IndexMat> tIds     = tInterp->get_ids();
                    // iterate through basis functions
                    for(moris::uint  iB = 0; iB <tOwners.numel(); iB++)
                    {
                        if(tOwners(iB) == par_rank())
                        {
                            tOwnedIndices.push_back(tIndices(iB));
                            tOwnedIds.push_back(tIds(iB));
                        }
                    }
                }
            }


            unique(tOwnedIndices);
            unique(tOwnedIds);


            std::cout<<"tOwnedIndices.size() = "<<tOwnedIndices.size()<<std::endl;
            std::cout<<"tOwnedIds.size() = "<<tOwnedIds.size()<<std::endl;

            delete tInterpMesh;

    }
}
}


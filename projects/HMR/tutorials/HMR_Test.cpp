/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Test.cpp
 *
 */

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Logger.hpp"
#include "cl_GlobalClock.hpp" // MRS/IOS/src
#include "cl_Tracer.hpp" // MRS/IOS/src

#include "cl_Matrix.hpp" //LINALG/src
#include "op_times.hpp" //LINALG/src
#include "fn_norm.hpp"

#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#define protected public
#define private   public
#include "cl_HMR_T_Matrix.hpp" //HMR/src
#undef protected
#undef private

// select namespaces
using namespace moris;
using namespace hmr;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;
//------------------------------------------------------------------------------

int
main(
        int    argc,
        char * argv[] )
{
//    // initialize MORIS global communication manager
//    gMorisComm = moris::Comm_Manager( &argc, &argv );
//
//    // Severity level 0 - all outputs
//    gLogger.initialize( 0 );
////------------------------------------------------------------------------------
//    // create settings object
//                moris::hmr::Parameters * tParameters = new moris::hmr::Parameters;
//
//                // this geometry creates one element, geometry coordinates
//                // are identical to parameter coordinates
//                moris::Matrix< moris::DDLUMat > tNumberOfElements = { {1}, {1}, {1} };
//                tParameters->set_number_of_elements_per_dimension( tNumberOfElements );
//
//                // parameter space goes from -1 to +1
//                moris::Matrix< moris::DDRMat > tDomain = { {2}, {2}, {2} };
//                tParameters->set_domain_dimensions( tDomain );
//                moris::Matrix< moris::DDRMat > tOffset = { {-1}, {-1}, {-1} };
//                tParameters->set_domain_offset( tOffset );
//
//                for( uint tOrder=1; tOrder<=3; ++tOrder )
//                {
//                    // set buffer size
//                    tParameters->set_refinement_buffer( tOrder );
//                    tParameters->set_staircase_buffer( tOrder );
//
//                    // activate truncation
//                    tParameters->set_bspline_truncation( true );
//
//                    // create factory
//                    Factory tFactory( mParameters );
//
//                    // create background mesh object
//                    moris::hmr::Background_Mesh_Base* tBackgroundMesh
//                    = tFactory.create_background_mesh( tParameters );
//
//
//                    // create B-Spline Mesh
//                    moris::hmr::BSpline_Mesh_Base* tBSplineMesh
//                    =  tFactory.create_bspline_mesh(
//                            tParameters,
//                            tBackgroundMesh,
//                            tParameters->get_bspline_input_pattern(),
//                            tOrder );
//
//                   moris::Cell< moris::hmr::BSpline_Mesh_Base*  > tBSplineMeshes( 1, tBSplineMesh );
//
//                    // create B-Spline Mesh
//                    moris::hmr::Lagrange_Mesh_Base* tLagrangeMesh
//                    =  tFactory.create_lagrange_mesh(
//                            tParameters,
//                            tBackgroundMesh,
//                            tBSplineMeshes,
//                            tParameters->get_lagrange_input_pattern(),
//                            tOrder );
//
//                    // create T-Matrix object
//                    moris::hmr::T_Matrix * tTMatrix
//                    = new moris::hmr::T_Matrix( tParameters, tLagrangeMesh, tBSplineMesh );
//
//                    // ask Lagrange mesh for number of nodes per element
//                    moris::luint tNumberOfNodes
//                    = tLagrangeMesh->get_number_of_nodes_on_proc();
//
//                    // this flag must stay true all the time
//                    bool tCheck = true;
//
//                    // loop over all nodes
//                    for( uint k=0; k<tNumberOfNodes; ++k )
//                    {
//                        // get pointer to node
//                        moris::hmr::Basis * tNode
//                        = tLagrangeMesh->get_basis_by_memory_index( k );
//
//                        // get node coordinate
//                        auto tXYZ = tNode->get_coords();
//
//                        // shape function vector
//                        moris::Matrix< moris::DDRMat > tN( tNumberOfNodes, 1 );
//
//                        tTMatrix->lagrange_shape_3d( tXYZ, tN );
//
//                        // epsilon environment
//                        moris::real tEpsilon = 1e-12;
//
//                        for( uint i=0; i<tNumberOfNodes; ++i )
//                        {
//                            if( i == k )
//                            {
//                                tCheck = tCheck && std::abs( tN( i ) - 1 ) < tEpsilon;
//                            }
//                            else
//                            {
//                                tCheck = tCheck && std::abs( tN( i ) )  < tEpsilon;
//                            }
//                        }
//                    }
//
//                    // tCheck must be true in order to pass test
//                    REQUIRE( tCheck );
//
//                    // tidy up memory
//                    delete tTMatrix;
//                    delete tBSplineMesh;
//                    delete tLagrangeMesh;
//                    delete tBackgroundMesh;
//                }
//
//                // delete settings object
//                delete tParameters;
//
//
////------------------------------------------------------------------------------
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
//    return 0;
}


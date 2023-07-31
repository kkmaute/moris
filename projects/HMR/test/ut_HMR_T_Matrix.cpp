/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_T_Matrix_Private.cpp
 *
 */

#include <catch.hpp>

#include "cl_HMR_T_Matrix.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src

#include "cl_Matrix.hpp" //LINALG/src
#include "op_times.hpp" //LINALG/src

namespace moris::hmr
{
    // Test class for T-matrix
    template< uint N >
    class T_Matrix_Test : public T_Matrix< N >
    {
    public:

        // Constructor
        using T_Matrix< N >::T_Matrix;

        // Test evaluation
        void evaluate_shape_function_test( const Matrix< DDRMat >& aXi, Matrix< DDRMat >& aN )
        {
            T_Matrix< N >::evaluate_shape_function( aXi, aN );
        }
    };

    // -----------------------------------------------------------------------------------------------------------------

    TEST_CASE( "HMR_T_Matrix_Private", "[moris],[mesh],[hmr],[hmr_t_matrix]" )
    {
        // these tests are only performed in serial. They have nothing to do with parallel.
        if( moris::par_size() == 1 )
        {

            // This test checks that the Lagrange shape functions work as expected
            SECTION ( "T-Matrix: Lagrange shape 2D" )
            {
                // create settings object
                auto tParameters = new moris::hmr::Parameters;

                // this geometry creates one element, geometry coordinates
                // are identical to parameter coordinates
                tParameters->set_number_of_elements_per_dimension( { {1}, {1} } );

                // parameter space goes from -1 to +1
                tParameters->set_domain_dimensions( { {2}, {2} } );
                tParameters->set_domain_offset( { {-1}, {-1} } );

                for( uint tOrder = 1; tOrder <= 3; ++tOrder )
                {
                    // set buffer size
                    tParameters->set_refinement_buffer( tOrder );
                    tParameters->set_staircase_buffer( tOrder );

                    // activate truncation
                    tParameters->set_bspline_truncation( true );

                    // create factory
                    moris::hmr::Factory tFactory( tParameters );

                    // create background mesh object
                    moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

                    // create B-Spline Mesh
                    moris::hmr::BSpline_Mesh_Base * tBSplineMesh =  tFactory.create_bspline_mesh(
                            tBackgroundMesh,
                            tParameters->get_bspline_input_pattern(),
                            tOrder );

                    // create container of B-Spline meshes
                    moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes( 1, tBSplineMesh );

                    // create B-Spline Mesh
                    moris::hmr::Lagrange_Mesh_Base * tLagrangeMesh =  tFactory.create_lagrange_mesh(
                            tBackgroundMesh,
                            tBSplineMeshes,
                            tParameters->get_lagrange_input_pattern(),
                            tOrder );

                    // create T-Matrix object
                    auto tTMatrix = new T_Matrix_Test< 2 >( tLagrangeMesh, tBSplineMesh );

                    // ask Lagrange mesh for number of nodes per element
                    moris::luint tNumberOfNodes = tLagrangeMesh->get_number_of_nodes_on_proc();

                    // this flag must stay true all the time
                    bool tCheck = true;

                    // shape function vector
                    moris::Matrix< moris::DDRMat > tN( tNumberOfNodes, 1 );

                    // loop over all nodes
                    for( uint k = 0; k < tNumberOfNodes; ++k )
                    {
                        // get pointer to node
                        moris::hmr::Basis * tNode = tLagrangeMesh->get_basis_by_memory_index( k );

                        // get node coordinate
                        auto tXY = tNode->get_coords();

                        tTMatrix->evaluate_shape_function_test( tXY, tN );

                        // epsilon environment
                        moris::real tEpsilon = 1e-12;

                        for( uint i = 0; i < tNumberOfNodes; ++i )
                        {
                            if( i == k )
                            {
                                tCheck = tCheck && std::abs( tN( i ) - 1 ) < tEpsilon;
                            }
                            else
                            {
                                tCheck = tCheck && std::abs( tN( i ) )  < tEpsilon;
                            }
                        }
                    }

                    // tCheck must be true in order to pass test
                    REQUIRE( tCheck );

                    // tidy up memory
                    delete tTMatrix;
                    delete tBSplineMesh;
                    delete tLagrangeMesh;
                    delete tBackgroundMesh;
                }

                // delete settings object
                delete tParameters;
            } //end section

    //-------------------------------------------------------------------------------

            // This test checks that the Lagrange shape functions work as expected
            SECTION ( "T-Matrix: Lagrange shape 3D" )
            {
                // create settings object
                auto tParameters = new moris::hmr::Parameters;

                // this geometry creates one element, geometry coordinates
                // are identical to parameter coordinates
                moris::Matrix< moris::DDLUMat > tNumberOfElements = { {1}, {1}, {1} };
                tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

                // parameter space goes from -1 to +1
                moris::Matrix< moris::DDRMat > tDomain = { {2}, {2}, {2} };
                tParameters->set_domain_dimensions( tDomain );
                moris::Matrix< moris::DDRMat > tOffset = { {-1}, {-1}, {-1} };
                tParameters->set_domain_offset( tOffset );

                for( uint tOrder=1; tOrder<=3; ++tOrder )
                {
                    // set buffer size
                    tParameters->set_refinement_buffer( tOrder );
                    tParameters->set_staircase_buffer( tOrder );

                    // activate truncation
                    tParameters->set_bspline_truncation( true );

                    // create factory
                    moris::hmr::Factory tFactory( tParameters );

                    // create background mesh object
                    moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

                    // create B-Spline Mesh
                    moris::hmr::BSpline_Mesh_Base * tBSplineMesh =  tFactory.create_bspline_mesh(
                            tBackgroundMesh,
                            tParameters->get_bspline_input_pattern(),
                            tOrder );

                   moris::Cell< moris::hmr::BSpline_Mesh_Base*  > tBSplineMeshes( 1, tBSplineMesh );

                    // create B-Spline Mesh
                    moris::hmr::Lagrange_Mesh_Base * tLagrangeMesh =  tFactory.create_lagrange_mesh(
                            tBackgroundMesh,
                            tBSplineMeshes,
                            tParameters->get_lagrange_input_pattern(),
                            tOrder );

                    // create T-Matrix object
                    auto tTMatrix = new T_Matrix_Test< 3 >( tLagrangeMesh, tBSplineMesh );

                    // ask Lagrange mesh for number of nodes per element
                    moris::luint tNumberOfNodes = tLagrangeMesh->get_number_of_nodes_on_proc();

                    // this flag must stay true all the time
                    bool tCheck = true;

                    // loop over all nodes
                    for( uint k=0; k<tNumberOfNodes; ++k )
                    {
                        // get pointer to node
                        moris::hmr::Basis * tNode = tLagrangeMesh->get_basis_by_memory_index( k );

                        // get node coordinate
                        auto tXYZ = tNode->get_coords();

                        // shape function vector
                        moris::Matrix< moris::DDRMat > tN( tNumberOfNodes, 1 );

                        tTMatrix->evaluate_shape_function_test( tXYZ, tN );

                        // epsilon environment
                        moris::real tEpsilon = 1e-12;

                        for( uint i=0; i<tNumberOfNodes; ++i )
                        {
                            if( i == k )
                            {
                                tCheck = tCheck && std::abs( tN( i ) - 1 ) < tEpsilon;
                            }
                            else
                            {
                                tCheck = tCheck && std::abs( tN( i ) )  < tEpsilon;
                            }
                        }
                    }

                    // tCheck must be true in order to pass test
                    REQUIRE( tCheck );

                    // tidy up memory
                    delete tTMatrix;
                    delete tBSplineMesh;
                    delete tLagrangeMesh;
                    delete tBackgroundMesh;
                }

                // delete settings object
                delete tParameters;
            } // end section
        } // end par rank
    } // end test case
}

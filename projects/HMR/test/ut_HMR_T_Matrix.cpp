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
        void evaluate_shape_function_test(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >& aN )
        {
            T_Matrix< N >::evaluate_shape_function( aXi, aN );
        }
    };

    // -----------------------------------------------------------------------------------------------------------------

    TEST_CASE( "HMR T-matrix", "[moris],[mesh],[hmr],[hmr_t_matrix]" )
    {
        // these tests are only performed in serial. They have nothing to do with parallel.
        if ( par_size() == 1 )
        {
            // Create parameters
            auto tParameters = new Parameters;

            // Loop over number of dimensions
            for ( uint iNumberOfDimensions = 2; iNumberOfDimensions <= 3; iNumberOfDimensions++ )
            {
                // Create one element
                tParameters->set_number_of_elements_per_dimension( Matrix< DDLUMat >( iNumberOfDimensions, 1, 1 ) );

                // Domain from -1 to 1 in each dimension
                tParameters->set_domain_dimensions( Matrix< DDRMat >( iNumberOfDimensions, 1, 2 ) );
                tParameters->set_domain_offset( Matrix< DDRMat >( iNumberOfDimensions, 1, -1 ) );

                // Loop over orders
                for ( uint iOrder = 1; iOrder <= 3; iOrder++ )
                {
                    // set buffer size
                    tParameters->set_refinement_buffer( iOrder );
                    tParameters->set_staircase_buffer( iOrder );

                    // activate truncation
                    tParameters->set_bspline_truncation( true );

                    // create factory
                    Factory tFactory( tParameters );

                    // create background mesh object
                    Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

                    // create B-Spline Mesh
                    BSpline_Mesh_Base* tBSplineMesh = tFactory.create_bspline_mesh(
                            tBackgroundMesh,
                            tParameters->get_bspline_input_pattern(),
                            iOrder );

                    // create container of B-Spline meshes
                    Cell< BSpline_Mesh_Base* > tBSplineMeshes( 1, tBSplineMesh );

                    // create B-Spline Mesh
                    Lagrange_Mesh_Base* tLagrangeMesh = tFactory.create_lagrange_mesh(
                            tBackgroundMesh,
                            tBSplineMeshes,
                            tParameters->get_lagrange_input_pattern(),
                            iOrder );

                    // create T-Matrix object
                    T_Matrix_Base* tTMatrix;
                    if ( iNumberOfDimensions == 2 )
                    {
                        tTMatrix = new T_Matrix_Test< 2 >( tLagrangeMesh, tBSplineMesh );
                    }
                    else
                    {
                        tTMatrix = new T_Matrix_Test< 3 >( tLagrangeMesh, tBSplineMesh );
                    }

                    // ask Lagrange mesh for number of nodes per element
                    luint tNumberOfNodes = tLagrangeMesh->get_number_of_nodes_on_proc();

                    // shape function vector
                    Matrix< DDRMat > tN( tNumberOfNodes, 1 );

                    // loop over all nodes
                    for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; iNodeIndex++ )
                    {
                        // get pointer to node
                        Basis* tNode = tLagrangeMesh->get_basis_by_memory_index( iNodeIndex );

                        // get node coordinate
                        auto tXYZ = tNode->get_coords();

                        // Evaluate shape function
                        if ( iNumberOfDimensions == 2 )
                        {
                            static_cast< T_Matrix_Test< 2 >* >( tTMatrix )->evaluate_shape_function_test( tXYZ, tN );
                        }
                        else
                        {
                            static_cast< T_Matrix_Test< 3 >* >( tTMatrix )->evaluate_shape_function_test( tXYZ, tN );
                        }

                        // epsilon environment
                        real tEpsilon = 1e-12;

                        // Check shape function vector
                        for ( uint iShapeFunctionIndex = 0; iShapeFunctionIndex < tNumberOfNodes; iShapeFunctionIndex++ )
                        {
                            if ( iShapeFunctionIndex == iNodeIndex )
                            {
                                CHECK( std::abs( tN( iShapeFunctionIndex ) - 1 ) < tEpsilon );
                            }
                            else
                            {
                                CHECK( std::abs( tN( iShapeFunctionIndex ) ) < tEpsilon );
                            }
                        }
                    }

                    // tidy up memory
                    delete tTMatrix;
                    delete tBSplineMesh;
                    delete tLagrangeMesh;
                    delete tBackgroundMesh;
                }
            }

            // delete settings object
            delete tParameters;
        }
    }
}

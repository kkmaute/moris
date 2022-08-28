/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix_2.cpp
 *
 */

#include "cl_HMR_T_Matrix_2.hpp" //HMR/src

#include <limits>

#include "HMR_Globals.hpp"     //HMR/src
#include "HMR_Tools.hpp"
#include "fn_eye.hpp"
#include "fn_norm.hpp"         //LINALG/src
#include "fn_sum.hpp"          //LINALG/src
#include "fn_trans.hpp"        //LINALG/src
#include "fn_inv.hpp"          //LINALG/src
#include "op_plus.hpp"         //LINALG/src
#include "op_times.hpp"        //LINALG/src

namespace moris
{
    namespace hmr
    {

        //-------------------------------------------------------------------------------

        T_Matrix_2::T_Matrix_2(
                const Parameters   * aParameters,
                BSpline_Mesh_Base  * aBSplineMesh,
                Lagrange_Mesh_Base * aLagrangeMesh,
                Lagrange_Mesh_Base * aLagrangeMeshCoarse)
        : T_Matrix( aParameters, aBSplineMesh, aLagrangeMesh ),
          mLagrangeMeshCoarse( aLagrangeMeshCoarse )
        {

        }

        //-------------------------------------------------------------------------------

        T_Matrix_2::~T_Matrix_2()
        {
        }

        //-------------------------------------------------------------------------------

        void T_Matrix_2::evaluate(
                const uint aBSplineMeshIndex,
                const bool aBool  )
        {
            // get B-Spline pattern of this mesh
            auto tBSplinePattern = mBSplineMesh->get_activation_pattern();

            // select pattern
            mLagrangeMesh->select_activation_pattern();

            // get number of elements on this Lagrange mesh
            luint tNumberOfElements = mLagrangeMesh->get_number_of_elements();

            // unflag all nodes on this mesh
            mLagrangeMeshCoarse->unflag_all_basis();
            uint tCoarseOrder = mLagrangeMeshCoarse->get_order();

            // number of nodes per element
            uint tNumberOfNodesPerElement = mLagrangeMesh->get_number_of_basis_per_element();
            uint tNumberOfCoarseNodesPerElement = mLagrangeMeshCoarse->get_number_of_basis_per_element();

            // unity matrix
            Matrix< DDRMat > tEye;
            eye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, tEye );

            // calculate transposed Lagrange T-Matrix
            Matrix< DDRMat > tL( this->get_lagrange_matrix() );

            /*if( mLagrangeMesh->get_activation_pattern()
                    == mParameters->get_lagrange_output_pattern() )
            {
                mLagrangeMesh->save_to_vtk("Lagrange.vtk");
                mBSplineMesh->save_to_vtk( "BSplines.vtk" );
                mLagrangeMesh->get_background_mesh()->save_to_vtk( "Background.vtk");
            }*/

            // loop over all elements
            for( luint e = 0; e < tNumberOfElements; ++e )
            {
                // get pointer to element
                auto tLagrangeElement       = mLagrangeMesh->get_element( e );
                auto tLagrangeElementCoarse = mLagrangeMeshCoarse->get_element( e );

                // FIXME : activate this flag
                //if ( tLagrangeElement->get_t_matrix_flag() )
                {
                    // get pointer to background element
                    auto tBackgroundElement = tLagrangeElement->get_background_element();

                    // initialize refinement Matrix
                    Matrix< DDRMat > tR;

                    bool tLagrangeEqualBspline = false;
                    bool tFirstLagrangeRefinement = true;

                    while( ! tBackgroundElement->is_active( tBSplinePattern ) )
                    {
                        if( tFirstLagrangeRefinement )
                        {
                            // right multiply refinement matrix
                            tR = this->get_refinement_matrix( tBackgroundElement->get_child_index() );
                            tFirstLagrangeRefinement = false;
                            tLagrangeEqualBspline = true;
                        }
                        else
                        {
                            tR = tR* this->get_refinement_matrix( tBackgroundElement->get_child_index() );
                        }

                        // jump to parent
                        tBackgroundElement = tBackgroundElement->get_parent();
                    }

                    // calculate the B-Spline T-Matrix
                    Matrix< DDRMat > tB;
                    Cell< Basis* > tDOFs;

                    this->calculate_t_matrix(
                            tBackgroundElement->get_memory_index(),
                            tB,
                            tDOFs );

                    // get change order matrix
                    const Matrix< DDRMat > & tC = this->get_change_order_matrix( tCoarseOrder );

                    Matrix< DDRMat > tT;

                    if(tLagrangeEqualBspline)
                    {
                        // transposed T-Matrix
                        tT = tC * tR * tL * tB;
                    }
                    else
                    {
                        tT = tC * tL * tB;
                    }

                    // number of columns in T-Matrix
                    uint tNCols = tT.n_cols();

                    // epsilon to count T-Matrix
                    real tEpsilon = 1e-12;

                    // loop over all nodes of this element
                    for( uint k = 0; k<tNumberOfCoarseNodesPerElement; ++k  )
                    {
                        // pointer to node
                        auto tNode = tLagrangeElementCoarse->get_basis( k );

                        // test if node is flagged
                        if ( ! tNode->is_flagged() )
                        {
                            // initialize counter
                            uint tCount = 0;

                            // reserve DOF cell
                            Cell< mtk::Vertex* > tNodeDOFs( tNCols, nullptr );

                            // reserve matrix with coefficients
                            Matrix< DDRMat > tCoefficients( tNCols, 1 );

                            // loop over all nonzero entries
                            for( uint i=0; i<tNCols; ++i )
                            {
                                if ( std::abs( tT( k, i ) ) > tEpsilon )
                                {
                                    // copy entry of T-Matrix
                                    tCoefficients( tCount ) = tT( k, i );

                                    // copy pointer of dof and convert to mtk::Vertex
                                    tNodeDOFs( tCount ) = tDOFs( i );

                                    // flag this DOF
                                    tDOFs( i )->flag();

                                    // increment counter
                                    ++tCount;
                                }
                            }

                            tCoefficients.resize( tCount, 1 );
                            tNodeDOFs.resize( tCount );

                            if ( aBool == true )
                            {
                                // init interpolation container for this node
                                tNode->init_interpolation( aBSplineMeshIndex );

                                // store the coefficients
                                tNode->set_weights( aBSplineMeshIndex, tCoefficients );

                                // store pointers to the DOFs
                                tNode->set_coefficients( aBSplineMeshIndex, tNodeDOFs );
                            }
                            // flag this node as processed
                            tNode->flag();
                        }
                    } // end loop over all nodes of this element
                } // end if element is flagged
            } // end loop over all elements
        }

        //-------------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */


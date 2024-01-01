/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix_Advanced.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_T_MATRIX_2_HPP_
#define SRC_HMR_CL_HMR_T_MATRIX_2_HPP_

#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "moris_typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "cl_Vector.hpp" //CNT/src
#include "cl_HMR_T_Matrix.hpp" //CNT/src

namespace moris::hmr
{
    /**
     * Advanced T-matrix class
     *
     * @tparam N Number of dimensions
     */
    template< uint N >
    class T_Matrix_Advanced : public T_Matrix< N >
    {
        //! Pointer to coarse Lagrange mesh
        Lagrange_Mesh_Base * mLagrangeMeshCoarse;

    public:

        /**
         * Constructor
         *
         * @param aLagrangeMesh Lagrange mesh pointer (fine)
         * @param aBSplineMesh B-spline Mesh pointer
         * @param aLagrangeMeshCoarse Lagrange mesh pointer (coarse)
         * @param aTruncate Whether or not to truncate B-splines
         */
        T_Matrix_Advanced(
                Lagrange_Mesh_Base * aLagrangeMesh,
                BSpline_Mesh_Base  * aBSplineMesh,
                Lagrange_Mesh_Base * aLagrangeMeshCoarse,
                bool                 aTruncate = true )
                : T_Matrix< N >( aLagrangeMesh, aBSplineMesh, aTruncate )
                , mLagrangeMeshCoarse( aLagrangeMeshCoarse )
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        ~T_Matrix_Advanced()
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        /**
         * Advanced evaluate function
         *
         * @param aBSplineMeshIndex B-spline mesh index
         * @param aBool FIXME what is this
         */
        void evaluate(
                uint aBSplineMeshIndex,
                bool aBool ) override
        {
            // get B-Spline pattern of this mesh
            auto tBSplinePattern = this->mBSplineMesh->get_activation_pattern();

            // select pattern
            this->mLagrangeMesh->select_activation_pattern();

            // get number of elements on this Lagrange mesh
            luint tNumberOfElements = this->mLagrangeMesh->get_number_of_elements();

            // unflag all nodes on this mesh
            this->mLagrangeMeshCoarse->unflag_all_basis();
            uint tCoarseOrder = this->mLagrangeMeshCoarse->get_order();

            // number of nodes per element
            uint tNumberOfNodesPerElement = this->mLagrangeMesh->get_number_of_bases_per_element();
            uint tNumberOfCoarseNodesPerElement = this->mLagrangeMeshCoarse->get_number_of_bases_per_element();

            // unity matrix
            Matrix< DDRMat > tEye;
            eye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, tEye );

            // calculate transposed Lagrange T-Matrix
            Matrix< DDRMat > tL( this->get_lagrange_matrix() );

            // loop over all elements
            for( luint e = 0; e < tNumberOfElements; ++e )
            {
                // get pointer to element
                auto tLagrangeElement       = this->mLagrangeMesh->get_element( e );
                auto tLagrangeElementCoarse = this->mLagrangeMeshCoarse->get_element( e );

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
                    Vector< Basis* > tDOFs;

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
                            Vector< mtk::Vertex* > tNodeDOFs( tNCols, nullptr );

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

                            if ( aBool )
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

        // -------------------------------------------------------------------------------------------------------------
    };
}

#endif /* SRC_HMR_CL_HMR_T_MATRIX_2_HPP_ */


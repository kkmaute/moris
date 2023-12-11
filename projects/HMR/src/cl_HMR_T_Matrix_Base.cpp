/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix_Base.cpp
 *
 */

#include "cl_HMR_T_Matrix.hpp"    //HMR/src
#include "HMR_Globals.hpp"        //HMR/src
#include "fn_eye.hpp"
#include "fn_inv.hpp"             //LINALG/src

namespace moris::hmr
{

    //-------------------------------------------------------------------------------

    T_Matrix_Base::T_Matrix_Base(
            Lagrange_Mesh_Base* aLagrangeMesh,
            BSpline_Mesh_Base*  aBSplineMesh,
            bool                aTruncate )
            : mLagrangeMesh( aLagrangeMesh )
            , mBSplineMesh( aBSplineMesh )
            , mTruncate( aTruncate )
    {
        this->init_lagrange_coefficients();
    }

    //-------------------------------------------------------------------------------

    T_Matrix_Base::~T_Matrix_Base()
    {
    }

    //-------------------------------------------------------------------------------

    void
    T_Matrix_Base::calculate_t_matrix(
            luint             aElementMemoryIndex,
            Matrix< DDRMat >& aTMatrixTransposed,
            Cell< Basis_Function* >&   aDOFs )
    {
        if ( mTruncate )
        {
            this->calculate_truncated_t_matrix( aElementMemoryIndex, aTMatrixTransposed, aDOFs );
        }
        else
        {
            this->calculate_untruncated_t_matrix( aElementMemoryIndex, aTMatrixTransposed, aDOFs );
        }
    }

    //-------------------------------------------------------------------------------

    void
    T_Matrix_Base::calculate_untruncated_t_matrix(
            luint             aElementMemoryIndex,
            Matrix< DDRMat >& aTMatrixTransposed,
            Cell< Basis_Function* >&   aDOFs )
    {
        aDOFs.clear();

        Element* tElement = mBSplineMesh->get_element_by_memory_index( aElementMemoryIndex );

        // get level of element
        uint tLevel = tElement->get_level();

        // get number of basis per element
        uint tNumberOfBasisPerElement = mBSplineMesh->get_number_of_bases_per_element();

        // help index for total number of basis
        uint tMaxNumberOfBasis = ( tLevel + 1 ) * tNumberOfBasisPerElement;

        // allocate max memory for matrix
        aTMatrixTransposed.set_size( tNumberOfBasisPerElement, tMaxNumberOfBasis, 0 );

        // allocate max memory for Basis
        aDOFs.resize( tMaxNumberOfBasis, nullptr );

        // write unity into level matrix
        Matrix< DDRMat > tT( mEye );

        // counter for basis
        uint tBasisCount = 0;

        // get pointer to parent
        Element* tParent = tElement;

        // loop over all levels and assemble transposed T-Matrix
        for ( uint iLevelIndex = 0; iLevelIndex <= tLevel; iLevelIndex++ )
        {
            // copy basis indices
            for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBasisPerElement; iBasisIndex++ )
            {
                // get pointer to basis
                Basis_Function* tBasisFunction = tParent->get_basis_function( iBasisIndex );

                // test if basis is active
                if ( tBasisFunction->is_active() )
                {
                    // copy columns into matrix
                    aTMatrixTransposed.set_column( tBasisCount, tT.get_column( iBasisIndex ) );

                    // copy pointer to basis into output array
                    aDOFs( tBasisCount++ ) = tBasisFunction;
                }
            }

            // left-multiply T-Matrix with child matrix
            tT = tT * mChildMatrices( tParent->get_background_element()->get_child_index() );

            // jump to next
            tParent = mBSplineMesh->get_parent_of_element( tParent );
        }

        // Shrink memory to exact size
        aDOFs.resize( tBasisCount, nullptr );
        aTMatrixTransposed.resize( tNumberOfBasisPerElement, tBasisCount );
    }

    //-------------------------------------------------------------------------------

    void
    T_Matrix_Base::calculate_truncated_t_matrix(
            luint             aElementMemoryIndex,
            Matrix< DDRMat >& aTMatrixTransposed,
            Cell< Basis_Function* >&   aDOFs )
    {
        // Clear adofs
        aDOFs.clear();

        // Get element from memory
        Element* tElement = mBSplineMesh->get_element_by_memory_index( aElementMemoryIndex );

        // get level of element
        uint tLevel = tElement->get_level();

        // get number of basis functions per element
        uint tNumberOfBFsPerElement = mBSplineMesh->get_number_of_bases_per_element();

        // help index for total number of basis
        uint tNumberOfBFsOnAllLevels = ( tLevel + 1 ) * tNumberOfBFsPerElement;

        // initialize counters
        uint tDOFCount   = 0;
        uint tNumBFsOnElement = 0;

        // initialize container of basis functions interpolating into the element
        moris::Cell< Basis_Function* > tBFsInterpolatingIntoElement( tNumberOfBFsOnAllLevels, nullptr );

        // size T-matrix
        aTMatrixTransposed.set_size( tNumberOfBFsPerElement, tNumberOfBFsOnAllLevels, 0.0 );
        aDOFs.resize( tNumberOfBFsOnAllLevels, nullptr );

        // copy basis functions on lowest level into output
        for ( uint iElemLocalBF = 0; iElemLocalBF < tNumberOfBFsPerElement; iElemLocalBF++ )
        {
            Basis_Function* tBasisFunction      = tElement->get_basis_function( iElemLocalBF );
            tBFsInterpolatingIntoElement( tNumBFsOnElement++ ) = tBasisFunction;

            if ( tBasisFunction->is_active() )
            {
                aTMatrixTransposed.set_column( tDOFCount, mEye.get_column( iElemLocalBF ) );
                aDOFs( tDOFCount++ ) = tBasisFunction;
            }
        }

        // if the current element is refined from level zero
        if ( tLevel > 0 )
        {
            // Get parent of element
            Element* tParentElement = mBSplineMesh->get_parent_of_element( tElement );

            // loop over all basis functions interpolating into the parent element
            for ( uint iParentElemLocalBF = 0; iParentElemLocalBF < tNumberOfBFsPerElement; iParentElemLocalBF++ )
            {
                // get a pointer to the current basis function of the parent element
                Basis_Function* tParentElemBF = tParentElement->get_basis_function( iParentElemLocalBF );

                // test if basis function is active
                if ( tParentElemBF->is_active() )
                {
                    // Get number of children of basis function
                    uint tNumberOfChildrenOfBasis = tParentElemBF->get_number_of_children();

                    // loop over all children of this basis function
                    for ( uint iChildNumber = 0; iChildNumber < tNumberOfChildrenOfBasis; iChildNumber++ )
                    {
                        // get pointer to child of basis function
                        Basis_Function* tChildBF = tParentElemBF->get_child( iChildNumber );

                        // test if child BF exists
                        if ( tChildBF != nullptr )
                        {
                            // test if child is not active
                            if ( !tChildBF->is_active() && !tChildBF->is_refined() )
                            {
                                // get memory index of child
                                luint tChildBasisFunctionMemoryIndex = tChildBF->get_memory_index();

                                // search for the child BF in the refined input element
                                for ( uint iBasisSearchIndex = 0; iBasisSearchIndex < tNumberOfBFsPerElement; iBasisSearchIndex++ )
                                {
                                    if ( tBFsInterpolatingIntoElement( iBasisSearchIndex )->get_memory_index() == tChildBasisFunctionMemoryIndex )
                                    {
                                        // fixme: this operation is supposed to work the same way for both Armadillo and Eigen.
#ifdef MORIS_USE_EIGEN
                                        auto tVal = aTMatrixTransposed.get_column( tDOFCount ).matrix_data() + mTruncationWeights( iChildNumber ) * mEye.get_column( iBasisSearchIndex ).matrix_data();
                                        aTMatrixTransposed.set_column( tDOFCount, tVal );
#else
                                        // compute the T-matrix weights in terms of the weights for the basis functions interpolating into the refined element
                                        auto tCol = aTMatrixTransposed.matrix_data().col( tDOFCount );
                                        tCol += mTruncationWeights( iChildNumber ) * mEye.matrix_data().col( iBasisSearchIndex );
#endif
                                        break;

                                    } // end if: child is found
                                } // end for: search for child in element
                            } // end if: child is neither active nor refined (i.e. it is disabled)
                        } // end if: child BF of the basis function interpolating into the parent element actually exists
                    } // end for: each child BF of the basis function interpolating into the parent element

                    // If there is are weights associated with the parent basis function, then add it to the list of DOFs for the input element
                    if ( aTMatrixTransposed.get_column( tDOFCount ).min() < -gEpsilon || aTMatrixTransposed.get_column( tDOFCount ).max() > gEpsilon )
                    {
                        // copy parent index
                        aDOFs( tDOFCount++ ) = tParentElemBF;
                    }
                    else // otherwise, don't and just have zero weights
                    {
                        // reset this column to zero
                        aTMatrixTransposed.set_column( tDOFCount, mZero.get_column( 0 ) );
                    }

                }    // end if: tBasisFunction is active

                // add the truncated parent basis functions to the list of BFs interpolating into the element
                tBFsInterpolatingIntoElement( tNumBFsOnElement++ ) = tParentElemBF;

            } // end for: each basis function interpolating into the parent element

        } // end if: the current element is refined from level zero

        // Shrink output to exact size
        aTMatrixTransposed.resize( tNumberOfBFsPerElement, tDOFCount );
        aDOFs.resize( tDOFCount, nullptr );

    }    // end function: hmr::T_Matrix_Base::calculate_truncated_t_matrix()

    //-------------------------------------------------------------------------------

    void
    T_Matrix_Base::init_lagrange_coefficients()
    {
        // number of Lagrange nodes per direction
        mLagrangeOrder = mLagrangeMesh->get_order();

        // Step 2: first, we initialize the parameter coordinates

        uint tNumberOfNodes = mLagrangeOrder + 1;

        // matrix containing parameter coordinates for points
        Matrix< DDRMat > tXi( tNumberOfNodes, 1, 0.0 );

        // step width
        real tDeltaXi = 2.0 / mLagrangeOrder;

        // first value
        tXi( 0 ) = -1.0;

        // intermediate values
        for ( uint iOrder = 1; iOrder < mLagrangeOrder; iOrder++ )
        {
            tXi( iOrder ) = tXi( iOrder - 1 ) + tDeltaXi;
        }

        // last value
        tXi( mLagrangeOrder ) = 1.0;

        // Step 3: we build a Vandermonde matrix
        Matrix< DDRMat > tVandermonde( tNumberOfNodes, tNumberOfNodes, 0.0 );

        for ( uint iRowIndex = 0; iRowIndex < tNumberOfNodes; iRowIndex++ )
        {
            for ( uint iColIndex = 0; iColIndex < tNumberOfNodes; iColIndex++ )
            {
                tVandermonde( iRowIndex, iColIndex ) = std::pow( tXi( iRowIndex ), tNumberOfNodes - iColIndex - 1 );
            }
        }

        // invert the Vandermonde matrix and store coefficients
        mLagrangeCoefficients = inv( tVandermonde );
    }

    //-------------------------------------------------------------------------------

    void
    T_Matrix_Base::evaluate(
            const uint aBSplineMeshIndex,
            const bool aBool )
    {
        // get B-Spline pattern of this mesh
        auto tBSplinePattern = mBSplineMesh->get_activation_pattern();

        // select pattern
        mLagrangeMesh->select_activation_pattern();

        // get number of elements on this Lagrange mesh
        luint tNumberOfElements = mLagrangeMesh->get_number_of_elements();

        // unflag all nodes on this mesh
        mLagrangeMesh->unflag_all_basis();

        // number of nodes per element
        uint tNumberOfNodesPerElement = mLagrangeMesh->get_number_of_bases_per_element();

        // unity matrix
        Matrix< DDRMat > tEye;
        eye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, tEye );

        // calculate transposed Lagrange T-Matrix
        Matrix< DDRMat > tL( this->get_lagrange_matrix() );

        // loop over all elements
        for ( luint iElementIndex = 0; iElementIndex < tNumberOfElements; iElementIndex++ )
        {
            // get pointer to element
            auto tLagrangeElement = mLagrangeMesh->get_element( iElementIndex );

            // FIXME : activate this flag
            // if ( tLagrangeElement->get_t_matrix_flag() )
            {
                // get pointer to background element
                auto tBackgroundElement = tLagrangeElement->get_background_element();

                // initialize refinement Matrix
                Matrix< DDRMat > tR;

                bool tLagrangeNotEqualBspline = false;
                bool tFirstLagrangeRefinement = true;

                while ( !tBackgroundElement->is_active( tBSplinePattern ) )
                {
                    if ( tFirstLagrangeRefinement )
                    {
                        // right multiply refinement matrix
                        tR                       = this->get_refinement_matrix( tBackgroundElement->get_child_index() );
                        tFirstLagrangeRefinement = false;
                        tLagrangeNotEqualBspline = true;
                    }
                    else
                    {
                        tR = tR * this->get_refinement_matrix( tBackgroundElement->get_child_index() );
                    }

                    // jump to parent
                    tBackgroundElement = tBackgroundElement->get_parent();
                }

                // calculate the B-Spline T-Matrix
                Matrix< DDRMat > tB;
                Cell< Basis_Function* >   tDOFs;

                uint tElemMemIndex = tBackgroundElement->get_memory_index();
                this->calculate_t_matrix(
                        tElemMemIndex,
                        tB,
                        tDOFs );

                Matrix< DDRMat > tT;

                if ( tLagrangeNotEqualBspline )    // need to consider refinement
                {
                    // transposed T-Matrix
                    tT = tR * tL * tB;
                }
                else    // no refinement
                {
                    tT = tL * tB;
                }

                // number of columns in T-Matrix
                uint tNCols = tT.n_cols();

                // epsilon to count T-Matrix
                real tEpsilon = 1e-12;

                // loop over all nodes of this element
                for ( uint iLagNode = 0; iLagNode < tNumberOfNodesPerElement; iLagNode++ )
                {
                    // pointer to node
                    auto tNode = tLagrangeElement->get_basis_function( iLagNode );

                    // test if node is flagged
                    if ( !tNode->is_flagged() )
                    {
                        // initialize counter
                        uint tNodeCount = 0;

                        // reserve DOF cell
                        Cell< mtk::Vertex* > tNodeDOFs( tNCols, nullptr );

                        // reserve matrix with coefficients
                        Matrix< DDRMat > tCoefficients( tNCols, 1 );

                        // loop over all nonzero entries
                        for ( uint iBspBF = 0; iBspBF < tNCols; ++iBspBF )
                        {
                            // get the T-matrix entry
                            real tTMatEntry = tT( iLagNode, iBspBF );

                            // ignore entries close to zero
                            if ( std::abs( tTMatEntry ) > tEpsilon )
                            {
                                // copy entry of T-Matrix
                                tCoefficients( tNodeCount ) = tTMatEntry;

                                // copy pointer of dof and convert to mtk::Vertex
                                tNodeDOFs( tNodeCount++ ) = tDOFs( iBspBF );

                                // flag this DOF
                                tDOFs( iBspBF )->flag();
                            }
                        }

                        tCoefficients.resize( tNodeCount, 1 );
                        tNodeDOFs.resize( tNodeCount );

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
                }    // end loop over all nodes of this element
            }        // end if element is flagged
        }            // end loop over all elements
    }

    //-------------------------------------------------------------------------------

    void
    T_Matrix_Base::evaluate_trivial(
            const uint aBSplineMeshIndex,
            const bool aBool )
    {
        // select pattern
        mLagrangeMesh->select_activation_pattern();

        // get number of elements on this Lagrange mesh
        luint tNumberOfElements = mLagrangeMesh->get_number_of_elements();

        // unflag all nodes on this mesh
        mLagrangeMesh->unflag_all_basis();

        // number of nodes per element
        uint tNumberOfNodesPerElement = mLagrangeMesh->get_number_of_bases_per_element();

        // loop over all elements
        for ( luint iElementIndex = 0; iElementIndex < tNumberOfElements; iElementIndex++ )
        {
            // get pointer to element
            auto tLagrangeElement = mLagrangeMesh->get_element( iElementIndex );

            // loop over all nodes of this element
            for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodesPerElement; iNodeIndex++ )
            {
                // pointer to node
                auto tNode = tLagrangeElement->get_basis_function( iNodeIndex );

                // test if node is flagged
                if ( !tNode->is_flagged() )
                {
                    // reserve matrix with coefficients
                    Matrix< DDRMat > tCoefficients( 1, 1, 1.0 );

                    // reserve DOF cell
                    Cell< mtk::Vertex* > tNodeDOFs( 1, tNode );

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
            }    // end loop over all nodes of this element
        }        // end loop over all elements
    }

    //-------------------------------------------------------------------------------

    void
    T_Matrix_Base::evaluate_extended_t_matrix(
            Element*                                    aBsplineElement,
            Element*                                    aLagrangeElement,
            moris::Cell< moris::Cell< mtk::Vertex* > >& aBsplineBasis,
            moris::Cell< Matrix< DDRMat > >&            aWeights )
    {
        Background_Element_Base* aBSpBackgroundElement = aBsplineElement->get_background_element();

        // evaluate
        this->init_modified_lagrange_parameter_coordinates( aLagrangeElement, aBsplineElement );
        this->recompute_lagrange_matrix();

        // get B-Spline pattern of this mesh
        uint tBSplinePattern = mBSplineMesh->get_activation_pattern();

        // select pattern
        mLagrangeMesh->select_activation_pattern();

        // unflag all nodes on this mesh
        mLagrangeMesh->unflag_all_basis();

        // number of nodes per element
        uint tNumberOfNodesPerElement = mLagrangeMesh->get_number_of_bases_per_element();

        // unity matrix
        Matrix< DDRMat > tEye;
        eye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, tEye );

        // calculate transposed Lagrange T-Matrix
        Matrix< DDRMat > tL( mTMatrixLagrangeModified );

        // get pointer to background element
        auto tBackgroundElement = aLagrangeElement->get_background_element();

        // initialize refinement Matrix
        Matrix< DDRMat > tR;

        bool tLagrangeEqualBspline    = false;
        bool tFirstLagrangeRefinement = true;

        while ( !tBackgroundElement->is_active( tBSplinePattern ) )
        {
            if ( tFirstLagrangeRefinement )
            {
                // right multiply refinement matrix
                tR                       = this->get_refinement_matrix( tBackgroundElement->get_child_index() );
                tFirstLagrangeRefinement = false;
                tLagrangeEqualBspline    = true;
            }
            else
            {
                tR = tR * this->get_refinement_matrix( tBackgroundElement->get_child_index() );
            }

            // jump to parent
            tBackgroundElement = tBackgroundElement->get_parent();
        }

        // calculate the B-Spline T-Matrix
        Matrix< DDRMat > tB;
        Cell< Basis_Function* >   tDOFs;

        this->calculate_t_matrix(
                aBSpBackgroundElement->get_memory_index(),
                tB,
                tDOFs );

        Matrix< DDRMat > tT;

        if ( tLagrangeEqualBspline )
        {
            // transposed T-Matrix
            tT = tR * tL * tB;
        }
        else
        {
            tT = tL * tB;
        }

        // number of columns in T-Matrix
        uint tNCols = tT.n_cols();

        // epsilon to count T-Matrix
        real tEpsilon = 1e-12;

        // resize the arrays
        aBsplineBasis.resize( tNumberOfNodesPerElement );
        aWeights.resize( tNumberOfNodesPerElement );

        // loop over all nodes of this element
        for ( uint iLagNode = 0; iLagNode < tNumberOfNodesPerElement; ++iLagNode )
        {
            // initialize counter
            uint tCount = 0;

            // reserve DOF cell
            Cell< mtk::Vertex* > tNodeDOFs( tNCols, nullptr );

            // reserve matrix with coefficients
            Matrix< DDRMat > tCoefficients( tNCols, 1 );

            // loop over all nonzero entries
            for ( uint iBspBF = 0; iBspBF < tNCols; ++iBspBF )
            {
                if ( std::abs( tT( iLagNode, iBspBF ) ) > tEpsilon )
                {
                    // copy entry of T-Matrix
                    tCoefficients( tCount ) = tT( iLagNode, iBspBF );

                    // copy pointer of dof and convert to mtk::Vertex
                    tNodeDOFs( tCount ) = tDOFs( iBspBF );

                    // flag this DOF
                    tDOFs( iBspBF )->flag();

                    // increment counter
                    ++tCount;
                }
            }

            tCoefficients.resize( tCount, 1 );
            tNodeDOFs.resize( tCount );

            aBsplineBasis( iLagNode ) = tNodeDOFs;
            aWeights( iLagNode )      = tCoefficients;
        }
    }

    //-------------------------------------------------------------------------------

}    // namespace moris::hmr

/*
 * cl_HMR_T_Matrix.cpp
 *
 *  Created on: Jun 23, 2018
 *      Author: messe
 */
#include <limits>

#include "op_times.hpp"        //LINALG/src
#include "fn_norm.hpp"         //LINALG/src
#include "fn_sum.hpp"          //LINALG/src
#include "fn_trans.hpp"        //LINALG/src
#include "fn_inv.hpp"          //LINALG/src
#include "op_plus.hpp"         //LINALG/src
#include "op_times.hpp"        //LINALG/src
#include "HMR_Globals.hpp"     //HMR/src
#include "cl_HMR_T_Matrix.hpp" //HMR/src


namespace moris
{
    namespace hmr
    {

//-------------------------------------------------------------------------------

        T_Matrix::T_Matrix(
                const Parameters   * aParameters,
                BSpline_Mesh_Base  * aBSplineMesh,
                Lagrange_Mesh_Base * aLagrangeMesh ) :
                            mParameters ( aParameters ),
                            mBSplineMesh( aBSplineMesh ),
                            mLagrangeMesh( aLagrangeMesh )
        {
            this->init_basis_index();
            this->init_unity_matrix();
            this->init_child_matrices();
            this->init_truncation_weights();
            this->init_lagrange_parameter_coordinates();
            this->init_lagrange_matrix();
            this->init_lagrange_coefficients();

            switch( mParameters->get_number_of_dimensions() )
            {
                case( 2 ) :
                {
                    mEvalNGeo   = & this->N_quad4;
                    mEvalN      = & T_Matrix :: lagrange_shape_2d;
                    mGetCorners = & this->get_child_corner_nodes_2d;

                    break;
                }
                case( 3 ) :
                {
                    mEvalNGeo   = & this->N_hex8;
                    mEvalN      = & T_Matrix :: lagrange_shape_3d;
                    mGetCorners = & this->get_child_corner_nodes_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown number of dimensions");
                    break;
                }
            }

            this->init_lagrange_refinement_matrices();
            //this->init_gauss_points();
            //this->init_mass_matrices();


            // set function pointer
            if ( aParameters->truncate_bsplines() )
            {
                mTMatrixFunction =
                        & T_Matrix::calculate_truncated_t_matrix;
            }
            else
            {
                mTMatrixFunction =
                        & T_Matrix::calculate_untruncated_t_matrix;
            }
        }

//-------------------------------------------------------------------------------

        T_Matrix::~T_Matrix()
        {

        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::calculate_t_matrix(
                const luint    & aMemoryIndex,
                Matrix< DDRMat >    & aTMatrixTransposed,
                Cell< Basis* > & aDOFs  )
        {
            ( this->*mTMatrixFunction )(
                    aMemoryIndex,
                    aTMatrixTransposed,
                    aDOFs  );
        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::calculate_untruncated_t_matrix(
                const luint    & aMemoryIndex,
                Matrix< DDRMat >    & aTMatrixTransposed,
                Cell< Basis* > & aDOFs  )
        {
            aDOFs.clear();

            Element* tElement = mBSplineMesh->get_element_by_memory_index( aMemoryIndex );

            // get level of element
            auto tLevel = tElement->get_level();

            // get number of basis per element
            auto tNumberOfBasisPerElement
                = mBSplineMesh->get_number_of_basis_per_element();

            // counter for basis
            uint tCount = 0;

            // get pointer to parent
            Element* tParent = tElement;

            // STEP 2: count number of active basis
            for( int l=tLevel; l>=0; --l )
            {
                // loop over all basis of this element and count active ones
                for( uint k=0; k<tNumberOfBasisPerElement; ++k )
                {
                    // test if basis is active
                    if( tParent->get_basis( k )->is_active() )
                    {
                        // increment counter
                        ++tCount;
                    }
                }

                // jump to parent
                tParent = mBSplineMesh->get_parent_of_element( tParent );
            }

            // STEP 3: initialize memory for T-Matrix and DOF indices

            // allocate memory for matrix
            aTMatrixTransposed.set_size( tNumberOfBasisPerElement, tCount, 0 );

            // allocate memory for Basis
            aDOFs.resize( tCount, nullptr );

            // STEP 4: calculate non-truncated matrix

            // write unity into level matrix
            Matrix< DDRMat > tT( mEye );

            // reset counter
            tCount = 0;

            // reset parent
            tParent = tElement;

            // loop over all levels and assemble transposed T-Matrix
            for( int l=tLevel; l>=0; --l )
            {
                // copy basis indices
                for( uint k=0; k<tNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis* tBasis = tParent->get_basis( k );

                    // test if basis is active
                    if ( tBasis->is_active() )
                    {
                        // copy columns into matrix
                        aTMatrixTransposed.cols( tCount, tCount ) = tT.cols( k, k );

                        // copy pointer to basis into output array
                        aDOFs( tCount++ ) = tBasis;
                    }

                }

                // left-multiply T-Matrix with child matrix

                // @fixme: the following operation causes a weird error in valgrind
                //         "Uninitialised value was created by a stack allocation"
                //
                tT = tT *mChild( tParent->get_background_element()->get_child_index() );
                //
                // therefore, the multiplication is done manually, which will cost
                // computation time

                /*auto tA = tT;
                auto tB = mChild( tParent->get_background_element()->get_child_index() );
                tT.fill( 0.0 );
                for( uint k=0; k<tNumberOfBasisPerElement; ++k )
                {
                    for( uint j=0; j<tNumberOfBasisPerElement; ++j )
                    {
                        for( uint i=0; i<tNumberOfBasisPerElement; ++i )
                        {
                            tT( i,k ) += tA( i,j ) * tB( j, k );
                        }
                    }
                }*/

                // jump to next
                tParent = mBSplineMesh->get_parent_of_element( tParent );

            }

        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::calculate_truncated_t_matrix(
                const luint    & aMemoryIndex,
                Matrix< DDRMat >    & aTMatrixTransposed,
                Cell< Basis* > & aDOFs  )
        {
            aDOFs.clear();

            Element* tElement = mBSplineMesh->get_element_by_memory_index( aMemoryIndex );

            // get level of element
            auto tLevel = tElement->get_level();

            // get number of basis per element
            auto tNumberOfBasisPerElement
                      = mBSplineMesh->get_number_of_basis_per_element();

            // help index for total number of basis
            uint tNumberOfBasis = ( tLevel+1 ) * tNumberOfBasisPerElement;

            // reserve memory for full T-Matrix
            Matrix< DDRMat > tTmatrixTransposed(
                    tNumberOfBasisPerElement,
                    tNumberOfBasis, 0 );

            // get pointer to parent
            Element* tParent = tElement;

            // initialize counter
            uint tCount = 0;

            // write unity into level matrix
            Matrix< DDRMat > tT( mEye );

            // help index for T-Matrix copying
            uint tEnd =  tNumberOfBasisPerElement - 1 ;

            // container for basis levels
            Matrix< DDLUMat > tBasisIndices( tNumberOfBasis, 1 );

            // create full matrix
            for( int l=tLevel; l>=0; --l )
            {
                // copy T-Matrix
                tTmatrixTransposed.cols( tCount, tCount+tEnd )
                        = tT.cols( 0, tEnd );

                // copy basis
                for( uint k=0; k<tNumberOfBasisPerElement; ++k )
                {
                    // copy memory index of basis
                    tBasisIndices( tCount++ )
                            =  tParent->get_basis( k )->get_memory_index();
                }

                // left-multiply T-Matrix with child matrix

                // @fixme: the following operation causes a weird error in valgrind
                //         "Uninitialised value was created by a stack allocation"
                //
                tT = tT *mChild( tParent->get_background_element()->get_child_index() );
                //
                // therefore, the multiplication is done manually, which will cost
                // computation time

                /* auto tA = tT;
                auto tB = mChild( tParent->get_background_element()->get_child_index() );
                tT.fill( 0.0 );
                for( uint k=0; k<tNumberOfBasisPerElement; ++k )
                {
                    for( uint j=0; j<tNumberOfBasisPerElement; ++j )
                    {
                        for( uint i=0; i<tNumberOfBasisPerElement; ++i )
                        {
                            tT( i,k ) += tA( i,j ) * tB( j, k );
                        }
                    }
                } */

                // jump to next parent
                tParent = mBSplineMesh->get_parent_of_element( tParent );

            }

            // reserve memory for truncated matrix
            Matrix< DDRMat > tTMatrixTruncatedTransposed(
                    tNumberOfBasisPerElement,
                    tNumberOfBasis, 0 );

            // reset counter
            tCount = 0;

            Matrix< DDLUMat > tDOFs( tNumberOfBasis, 1 );

            // copy basis on lowest level into output
            for( uint k=0; k<tNumberOfBasisPerElement; ++k )
            {
                if (  mBSplineMesh
                        ->get_basis_by_memory_index( tBasisIndices( k ) )->is_active() )
                {
                    tTMatrixTruncatedTransposed.cols( tCount, tCount )
                        = tTmatrixTransposed.cols( k, k );

                    tDOFs( tCount++ ) = tBasisIndices( k );
                }
            }

            uint tK0 = 0;
            uint tK1 = tNumberOfBasisPerElement;
            uint tK2 = 2* tNumberOfBasisPerElement;

            // ask B-Spline mesh for number of children per basis
            auto tNumberOfChildrenPerBasis = mBSplineMesh->get_number_of_children_per_basis();

            // loop over higher levels
            for( int l=tLevel-1; l>=0; --l )
            {
                // loop over all basis of this level
                for( uint k=tK1; k<tK2; ++k )
                {
                    // get pointer to basis
                    Basis * tBasis = mBSplineMesh
                            ->get_basis_by_memory_index( tBasisIndices( k ) );

                    // test if basis is active
                    if ( tBasis->is_active() )
                    {
                        //tTMatrixTruncatedTransposed.cols( tCount, tCount )
                        //        = tTmatrixTransposed.cols( k, k );

                        // loop over all children of this basis
                        for( uint j=0; j<tNumberOfChildrenPerBasis; ++j )
                        {
                            // get pointer to child of basis
                            Basis * tChild = tBasis->get_child( j );

                            // test if child exists
                            if ( tChild != NULL )
                            {
                                // test if child is active
                                if ( !tChild->is_active( ) && !tChild->is_refined() )
                                //if ( tChild->is_active()  )
                                {
                                    // get memory index of child
                                    luint tIndex = tChild->get_memory_index();

                                    // search for child in element
                                    for( uint i=tK0; i<tK1; ++i )
                                    {
                                        if ( tBasisIndices( i ) == tIndex )
                                        {
                                            // subtract column from matrix
                                            tTMatrixTruncatedTransposed.cols( tCount, tCount )
                                                += mTruncationWeights( j )
                                                  * tTmatrixTransposed.cols( i, i );
                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        // copy parent index
                        tDOFs( tCount++ ) = tBasisIndices( k );
                    }
                }

                tK0 = tK1;
                tK1 = tK2;
                tK2 += tNumberOfBasisPerElement;
            }

            // remember new number of basis, we can drop the old one
            tNumberOfBasis = tCount;

            // test vector
            Matrix< DDRMat > tCol( tNumberOfBasisPerElement, 1 );

            // reset counter
            tCount = 0;

            // flags ( to avoid calculating the norm twice )
            Matrix< DDUMat > tUseColumn( tNumberOfBasis, 1, 0 );

            // count number of relevant entries
            for( uint k=0; k<tNumberOfBasis; ++k )
            {
                // copy column
                tCol.cols( 0, 0 ) = tTMatrixTruncatedTransposed.cols( k, k );

                // test if matrix is relevant
                if ( norm(tCol) > gEpsilon )
                {
                    tUseColumn( k ) = 1;
                    ++tCount;
                }
            }

            // assign memory for output matrix
            aTMatrixTransposed.set_size( tNumberOfBasisPerElement, tCount );

            aDOFs.resize( tCount, nullptr );

            // reset counter
            tCount = 0;

            // loop over matrix and copy relevant entries
            for( uint k=0; k<tNumberOfBasis; ++k )
            {
                // test if matrix is relevant
                if ( tUseColumn( k ) == 1 )
                {
                    aTMatrixTransposed.cols( tCount, tCount )
                            = tTMatrixTruncatedTransposed.cols( k, k );

                    // get pointer to basis from background mesh
                    aDOFs( tCount++ ) =
                            mBSplineMesh->get_basis_by_memory_index( tDOFs( k ) );
                }

            }

        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::init_basis_index()
        {


            // get dimensions
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            Background_Element_Base* tBackElement = this->create_background_element();

            // create a prototype of a B-Spline Element
            Element* tElement = mBSplineMesh->create_element( tBackElement );

            // B-Spline order
            mBSplineOrder = mBSplineMesh->get_order();

            // number of nodes per direction
            uint tNodesPerDirection = mBSplineOrder + 1;

            // calculate number of basis per element
            uint tNumberOfBasis
                = std::pow( tNodesPerDirection , tNumberOfDimensions );

            // initialize index matrix
            mBasisIndex.set_size( tNumberOfBasis, 1 );

            mBSplineIJK.set_size( tNumberOfDimensions, tNumberOfBasis );

            // loop over all basis
            if ( tNumberOfDimensions == 2)
            {
                // container for ijk position of basis
                luint tIJ[ 2 ];
                for( uint k=0; k<tNumberOfBasis; ++k )
                {
                    // get position from element
                    tElement->get_ijk_of_basis( k, tIJ );

                    // calculate index in matrix
                     uint tIndex = tIJ[ 0 ] + tIJ[ 1 ]*tNodesPerDirection;

                     mBasisIndex( tIndex ) = k;
                     mBSplineIJK( 0, k ) = tIJ[ 0 ];
                     mBSplineIJK( 1, k ) = tIJ[ 1 ];

                }
            }
            else if ( tNumberOfDimensions == 3)
            {
                // container for ijk position of basis
                luint tIJK[ 3 ];
                for( uint k=0; k<tNumberOfBasis; ++k )
                {
                    // get position from element
                    tElement->get_ijk_of_basis( k, tIJK );

                    // calculate index in matrix
                    uint tIndex =  tIJK[ 0 ]
                                 + tNodesPerDirection *
                                 ( tIJK[ 1 ] + tIJK[ 2 ]*tNodesPerDirection );

                    mBasisIndex( tIndex ) = k;

                    mBSplineIJK( 0, k ) = tIJK[ 0 ];
                    mBSplineIJK( 1, k ) = tIJK[ 1 ];
                    mBSplineIJK( 2, k ) = tIJK[ 2 ];
                }
            }

            // tidy up
            delete tElement;
            delete tBackElement;
        }

//-------------------------------------------------------------------------------

        Background_Element_Base*
        T_Matrix::create_background_element()
        {

            Background_Element_Base* aBackElement = nullptr;

            // create a prototype for a background element
            switch ( mParameters->get_number_of_dimensions() )
            {
                case( 2 ) :
                {
                    luint tIJ[ 2 ] = { 0, 0 };
                    aBackElement = new Background_Element< 2, 4, 8 >(
                            ( Background_Element_Base* ) nullptr,
                            0,
                            tIJ,
                            0 ,
                            ( uint ) 0,
                            ( uint ) 0,
                            ( uint ) gNoProcOwner );
                    break;
                }
                case( 3 ) :
                {
                    luint tIJK[ 3 ] = { 0, 0, 0 };
                    aBackElement = new Background_Element< 3, 8, 26 >(
                            ( Background_Element_Base* ) nullptr,
                            0,
                            tIJK,
                            0 ,
                            ( uint ) 0,
                            ( uint ) 0,
                            ( uint ) gNoProcOwner );
                    break;
                }
                default :
                {
                     MORIS_ERROR( false, "unknown number of dimensions." );
                }
            }

            return aBackElement;
        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::init_unity_matrix()
        {
            // get number of basis per element
            uint tNumberOfBasis
                 = std::pow( mBSplineMesh->get_order()+1,
                    mParameters->get_number_of_dimensions() );

            mEye.set_size( tNumberOfBasis, tNumberOfBasis, 0.0 );

            // loop over all positions
            for( uint i=0; i<tNumberOfBasis; ++i )
            {
                mEye( i, i ) = 1.0;
            }
        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::init_child_matrices()
        {
            // get order of mesh
            uint tOrder = mBSplineMesh->get_order();

            // create temporary matrix
            Matrix< DDRMat > tFactors( tOrder+1, tOrder+2, 0.0 );

            // number of nodes per direction
            uint n = tOrder+1;

            // weight factor
            real tWeight = 1.0/ std::pow( 2, tOrder );

            for( uint j=0; j<=n; ++j )
            {
                for( uint i=0; i<=tOrder; ++i )
                {
                    uint k = tOrder - 2*i + j ;
                    if ( k <= n )
                    {
                        tFactors( i, j ) =  tWeight*this->nchoosek( n, k );
                    }
                }
            }

            // left matrix
            Matrix< DDRMat > TL( n, n, 0.0 );
            TL.cols( 0, tOrder ) = tFactors.cols( 0, tOrder );

            // right matrix
            Matrix< DDRMat > TR( n, n, 0.0 );
            TR.cols( 0, tOrder ) = tFactors.cols( 1, n );

            // get number of dimensions from settings
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // determine number of children
            uint tNumberOfChildren = std::pow( 2, tNumberOfDimensions );

            // determine number of basis per element
            uint tNumberOfBasis = std::pow( tOrder+1, tNumberOfDimensions );

            // empty matrix
            Matrix< DDRMat > tEmpty ( tNumberOfBasis, tNumberOfBasis, 0.0 );

            // container for child relation matrices ( transposed! )
            mChild.resize( tNumberOfChildren, tEmpty );

            // populate child matrices
            if ( tNumberOfDimensions == 2 )
            {
                uint b=0;
                for( uint l=0; l<n; ++l )
                {
                    for( uint k=0; k<n; ++k )
                    {
                        uint a=0;
                        for( uint j=0; j<n; ++j )
                        {
                            for( uint i=0; i<n; ++i )
                            {
                                mChild( 0 )( mBasisIndex( a ), mBasisIndex( b ) )
                                        = TL( k,i )*TL( l, j );
                                mChild( 1 )( mBasisIndex( a ), mBasisIndex( b ) )
                                        = TR( k,i )*TL( l, j );
                                mChild( 2 )( mBasisIndex( a ), mBasisIndex( b ) )
                                        = TL( k,i )*TR( l, j );
                                mChild( 3 )( mBasisIndex( a ), mBasisIndex( b ) )
                                        = TR( k,i )*TR( l, j );
                                ++a;
                            }
                        }
                        ++b;
                    }
                }
            }
            else if ( tNumberOfDimensions == 3 )
            {
                uint b=0;
                for( uint p=0; p<n; ++p )
                {
                    for( uint q=0; q<n; ++q )
                    {
                        for( uint l=0; l<n; ++l )
                        {
                            uint a=0;
                            for( uint k=0; k<n; ++k )
                            {
                                for( uint j=0; j<n; ++j )
                                {
                                    for( uint i=0; i<n; ++i )
                                    {
                                        mChild( 0 )( mBasisIndex( a ), mBasisIndex( b ) )
                                                   = TL( l, i ) * TL( q, j ) * TL( p, k );
                                        mChild( 1 )( mBasisIndex( a ), mBasisIndex( b ) )
                                                   = TR( l, i ) * TL( q, j ) * TL( p, k );
                                        mChild( 2 )( mBasisIndex( a ), mBasisIndex( b ) )
                                                   = TL( l, i ) * TR( q, j ) * TL( p, k );
                                        mChild( 3 )( mBasisIndex( a ), mBasisIndex( b ) )
                                                   = TR( l, i ) * TR( q, j ) * TL( p, k );
                                        mChild( 4 )( mBasisIndex( a ), mBasisIndex( b ) )
                                                   = TL( l, i ) * TL( q, j ) * TR( p, k );
                                        mChild( 5 )( mBasisIndex( a ), mBasisIndex( b ) )
                                                   = TR( l, i ) * TL( q, j ) * TR( p, k );
                                        mChild( 6 )( mBasisIndex( a ), mBasisIndex( b ) )
                                                   = TL( l, i ) * TR( q, j ) * TR( p, k );
                                        mChild( 7 )( mBasisIndex( a ), mBasisIndex( b ) )
                                                   = TR( l, i ) * TR( q, j ) * TR( p, k );
                                        ++a;
                                    }
                                }
                            }
                            ++b;
                        }
                    }
                }
            }
        }
//-------------------------------------------------------------------------------

        void
        T_Matrix::init_truncation_weights()
        {
            // get order of mesh
            uint tOrder = mBSplineMesh->get_order();

            // number of children per direcion
            uint tNumberOfChildren = tOrder + 2;

            // matrix containing 1D weights
            Matrix< DDRMat > tWeights( tNumberOfChildren, 1 );

            // scale factor for 1D weights
            real tScale = 1.0/( (real) std::pow( 2, tOrder ) );

            // calculate 1D weights
            for( uint k=0; k<tNumberOfChildren; ++k )
            {
                tWeights( k ) = tScale*this->nchoosek( tOrder+1, k );
            }

            // get number of dimensions from settings
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // allocate weights
            mTruncationWeights.set_size(
                    std::pow( tNumberOfChildren, tNumberOfDimensions), 1 );

            if ( tNumberOfDimensions == 2 )
            {
                // init counter
                uint tCount = 0;

                // loop over all positions
                for( uint j=0; j<tNumberOfChildren; ++j )
                {
                    for( uint i=0; i<tNumberOfChildren; ++i )
                    {
                        mTruncationWeights( tCount++ )
                                = tWeights( i ) * tWeights( j );
                    }
                }
            }
            else if ( tNumberOfDimensions == 3 )
            {
                // init counter
                uint tCount = 0;

                // loop over all positions
                for( uint k=0; k<tNumberOfChildren; ++k )
                {
                    for( uint j=0; j<tNumberOfChildren; ++j )
                    {
                        for( uint i=0; i<tNumberOfChildren; ++i )
                        {
                            mTruncationWeights( tCount++ )
                                    = tWeights( i ) * tWeights( j ) * tWeights( k );
                        }
                    }
                }
            }
        }

//-------------------------------------------------------------------------------
        void
        T_Matrix::init_lagrange_parameter_coordinates()
        {

            // get dimensions from settings
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // create background element for reference
            Background_Element_Base* tBackElement = this->create_background_element();

            // number of nodes per direction
            uint tNodesPerDirection = mLagrangeMesh->get_order() + 1;

            // number of nodes
            mNumberOfNodes = std::pow( tNodesPerDirection, tNumberOfDimensions );

            // create a Lagrange element
            Element *tElement = mLagrangeMesh->create_element( tBackElement );

            // assign memory for parameter coordinates
            mLagrangeParam.set_size( tNumberOfDimensions, mNumberOfNodes );


            // scaling factor
            real tScale = 1.0/( ( real ) mLagrangeMesh->get_order() );

            // container for ijk position of basis
            luint tIJK[ 3 ];

            // ijk positions for reference Lagrange element
            mLagrangeIJK.set_size(  tNumberOfDimensions, mNumberOfNodes );

            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                // get position from element
                tElement->get_ijk_of_basis( k, tIJK );

                // save coordinate into memory
                for( uint i=0; i<tNumberOfDimensions; ++i )
                {
                    mLagrangeParam( i, k ) = 2*tScale*tIJK[ i ] - 1.0;
                    mLagrangeIJK( i, k )   = tIJK[ i ];
                }
            }

            // tidy up
            delete tElement;
            delete tBackElement;
        }

//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
        void
        T_Matrix::init_lagrange_matrix()
        {
            // get number of dimensions from settings
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // get number of basis per element of B-Spline mesh
            uint tNumberOfBasis = mBSplineIJK.n_cols();

            // get number of Lagrange nodes
            uint tNumberOfNodes = mLagrangeParam.n_cols();

            // initialize T-Matrix
            mTMatrixLagrange.set_size( tNumberOfNodes, tNumberOfBasis, 1 );

            // get order of B-Spline mesh
            uint tOrder = mBSplineMesh->get_order();

            // loop over all Lagrange nodes
            for( uint k = 0; k<tNumberOfNodes; ++k)
            {
                // loop over all B-Spline Basis
                for( uint j=0; j<tNumberOfBasis; ++j )
                {

                    // loop over all dimensions
                    for( uint i=0; i<tNumberOfDimensions; ++i )
                    {
                        mTMatrixLagrange( k, j ) *= this->b_spline_shape_1d(
                                tOrder,
                                mBSplineIJK( i, j ),
                                mLagrangeParam( i, k ) );
                    }
                }
            }
        }

//-------------------------------------------------------------------------------

        real
        T_Matrix::nchoosek( const uint & aN, const uint aK )
        {
            real aResult = 1.0;

            for ( uint i=1; i<=aK; ++i )
            {
                aResult *= ( ( real ) aN+1-i ) / ( real( i ) );
            }

            return aResult;
        }

//-------------------------------------------------------------------------------

        /**
         * 1D shape function
         */
        real
        T_Matrix::b_spline_shape_1d(
                const uint & aOrder,
                const uint & aK,
                const real & aXi ) const
        {
            // max number of entries in lookup table
            uint tSteps = 2*(aOrder + 1 );

            // temporary matrix that contains B-Spline segments
            Matrix< DDRMat > tDeltaXi( tSteps, 1, 0 );
            for( uint i=0; i< tSteps; ++i )
            {
                tDeltaXi( i ) = ( ( ( real ) i ) - ( ( real ) aOrder ) )*2.0 - 1.0;
            }

            // temporary matrix that contains evaluated values
            Matrix< DDRMat > tN( aOrder+1, 1, 0 );

            // initialize zero order values
            for( uint i=0; i<= aOrder; ++i )
            {
                if ( tDeltaXi( i + aK ) <= aXi && aXi < tDeltaXi( i + aK + 1 ) )
                {
                    tN( i ) = 1.0;
                }
            }

            // loop over all orders
            for( uint r=1; r<=aOrder; ++r )
            {
                // copy values of tN into old matrix
                Matrix< DDRMat > tNold( tN );

                // loop over all contributions
                for( uint i=0; i<=aOrder-r; ++i )
                {

                    // help values
                    real tA = aXi - tDeltaXi( i+aK );
                    real tB = tDeltaXi( i+aK + r + 1 ) - aXi;

                    tN( i ) = 0.5*( tA * tNold( i ) + tB * ( tNold( i + 1 ) ) ) / ( ( real ) r );
                }
            }

            // first value in entry is shape value
            return tN( 0 );

        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::b_spline_shape(
                const real        & aXi,
                const real        & aEta,
                Matrix< DDRMat >       & aN ) const
        {
            // evaluate contributions for xi and eta
            Matrix< DDRMat >  tNxi( mBSplineOrder+1, 1 );
            Matrix< DDRMat > tNeta( mBSplineOrder+1, 1 );

            for( uint i=0; i<=mBSplineOrder; ++i )
            {
                tNxi( i ) = this->b_spline_shape_1d(
                                mBSplineOrder,
                                i,
                                aXi );
            }
            for( uint j=0; j<=mBSplineOrder; ++j )
            {
                tNeta( j )  = this->b_spline_shape_1d(
                        mBSplineOrder,
                        j,
                        aEta );
            }

            // create shape vector in correct order
            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                aN( k ) =     tNxi( mBSplineIJK( 0, k ) )
                           * tNeta( mBSplineIJK( 1, k ) );
            }

        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::b_spline_shape(
                const real        & aXi,
                const real        & aEta,
                const real        & aZeta,
                Matrix< DDRMat >       & aN ) const
        {
            // evaluate contributions for xi and eta
            Matrix< DDRMat >  tNxi( mBSplineOrder+1, 1 );
            Matrix< DDRMat > tNeta( mBSplineOrder+1, 1 );
            Matrix< DDRMat > tNzeta( mBSplineOrder+1, 1 );

            for( uint i=0; i<=mBSplineOrder; ++i )
            {
                tNxi( i ) = this->b_spline_shape_1d(
                        mBSplineOrder,
                        i,
                        aXi );
            }

            for( uint j=0; j<=mBSplineOrder; ++j )
            {
                tNeta( j )  = this->b_spline_shape_1d(
                        mBSplineOrder,
                        j,
                        aEta );
            }

            for( uint k=0; k<=mBSplineOrder; ++k )
            {
                tNzeta( k )  = this->b_spline_shape_1d(
                        mBSplineOrder,
                        k,
                        aZeta );
            }


            // create shape vector in correct order
            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                aN( k ) =        tNxi( mBSplineIJK( 0, k ) )
                             *  tNeta( mBSplineIJK( 1, k ) )
                             * tNzeta( mBSplineIJK( 2, k ) );
            }

        }
//-------------------------------------------------------------------------------

        void
        T_Matrix::init_lagrange_coefficients()
        {
            // number of Lagrange nodes per direction
            mLagrangeOrder = mLagrangeMesh->get_order();

            // Step 2: first, we initialize the parameter coordinates

            uint tNumberOfNodes = mLagrangeOrder + 1;

            // matrix containing parameter coordinates for points
            Matrix< DDRMat > tXi( tNumberOfNodes, 1, 0.0 );

            // stepwidth
            real tDeltaXi = 2.0/mLagrangeOrder;

            // first value
            tXi( 0 ) = -1.0;

            // intermediate values
            for( uint k=1; k<mLagrangeOrder; ++k )
            {
                tXi( k ) = tXi( k-1 ) + tDeltaXi;
            }

            // last value
            tXi( mLagrangeOrder ) = 1.0;

            // Step 3: we build a Vandermonde matrix
            Matrix< DDRMat > tVandermonde( tNumberOfNodes, tNumberOfNodes, 0.0 );

            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                for( uint i=0; i<tNumberOfNodes; ++i )
                {
                    tVandermonde( k, i ) = std::pow(
                            tXi( k ),
                            tNumberOfNodes-i-1 );
                }
            }



            // invert the Vandermonde matrix and store coefficients
            mLagrangeCoefficients.set_size( tNumberOfNodes, tNumberOfNodes, 1 );

            Matrix< DDRMat > tIV = inv( tVandermonde );

            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                // help vector
                Matrix< DDRMat > tRHS( tNumberOfNodes, 1, 0.0 );
                tRHS( k ) = 1.0;
                Matrix< DDRMat > tLHS = tIV * tRHS;
                mLagrangeCoefficients.cols( k, k ) = tLHS.cols( 0, 0 );
            }
        }

//-------------------------------------------------------------------------------

        real
        T_Matrix::lagrange_shape_1d(
                const uint        & aBasisNumber,
                const real        & aXi ) const
        {
            // use horner scheme to evaluate 1D Lagrange function
            real aResult = 0.0;
            for( uint i=0; i<mLagrangeOrder; ++i )
            {
                aResult =
                        ( aResult + mLagrangeCoefficients( i, aBasisNumber ) )
                        *aXi;
            }

            return aResult + mLagrangeCoefficients( mLagrangeOrder, aBasisNumber );
        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::lagrange_shape_2d(
                const Matrix< DDRMat > & aXi,
                      Matrix< DDRMat > & aN ) const
        {
            // evaluate contributions for xi and eta
            Matrix< DDRMat >  tNxi( mLagrangeOrder+1, 1 );
            Matrix< DDRMat > tNeta( mLagrangeOrder+1, 1 );

            for( uint i=0; i<=mLagrangeOrder; ++i )
            {
                tNxi( i ) = this->lagrange_shape_1d( i, aXi( 0 ) );
            }
            for( uint j=0; j<=mLagrangeOrder; ++j )
            {
                tNeta( j ) = this->lagrange_shape_1d( j, aXi( 1 ) );
            }

            // create shape vector in correct order
            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                    aN( k ) =    tNxi( mLagrangeIJK( 0, k ) )
                              * tNeta( mLagrangeIJK( 1, k ) );
            }

        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::lagrange_shape_3d(
                const Matrix< DDRMat > & aXi,
                      Matrix< DDRMat > & aN ) const
        {
            // evaluate contributions for xi and eta and zeta
            Matrix< DDRMat >   tNxi( mLagrangeOrder+1, 1 );
            Matrix< DDRMat >  tNeta( mLagrangeOrder+1, 1 );
            Matrix< DDRMat > tNzeta( mLagrangeOrder+1, 1 );

            for( uint i=0; i<=mLagrangeOrder; ++i )
            {
                tNxi( i ) = this->lagrange_shape_1d( i, aXi( 0 ) );
            }

            for( uint j=0; j<=mLagrangeOrder; ++j )
            {
                tNeta( j ) = this->lagrange_shape_1d( j, aXi( 1 ) );
            }

            for( uint k=0; k<=mLagrangeOrder; ++k )
            {
                tNzeta( k ) = this->lagrange_shape_1d( k, aXi( 2 ) );
            }

            // create shape vector in correct order
            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                    aN( k ) =     tNxi( mLagrangeIJK( 0, k ) )
                              *  tNeta( mLagrangeIJK( 1, k ) )
                              * tNzeta( mLagrangeIJK( 2, k ) );
            }
        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::evaluate()
        {
            // get B-Spline pattern of this mesh
            auto tBSplinePattern = mBSplineMesh->get_activation_pattern();

            // select pattern
            mLagrangeMesh->select_activation_pattern();

            // get number of elements on this Lagrange mesh
            auto tNumberOfElements = mLagrangeMesh->get_number_of_elements();

            // unflag all bsplines
            //mBSplineMesh->unflag_all_basis();

            // flag B-Splines
            //uint tNumberOfBSplinesPerElement
            //    = mBSplineMesh->get_number_of_basis_per_element();

            // initialize counter
            /*luint tCount = 0;



            // loop over all elements
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to B-Spline Element
                auto tElement = mBSplineMesh->get_element( e );

                // check if Element is flagged
                if ( tElement->get_t_matrix_flag() )
                {
                    for ( uint k=0; k<tNumberOfBSplinesPerElement; ++k )
                    {

                        // get basis
                        auto tBasis = tElement->get_basis( k );

                        if ( ! tBasis->is_flagged() && tBasis->is_active() )
                        {
                            // set index for basis
                            tBasis->set_local_index( tCount++ );

                            // flag basis
                            tBasis->flag();
                        }

                    }
                }
                } */


            // unflag all nodes on this mesh
            mLagrangeMesh->unflag_all_basis();

            // number of nodes per element
            auto tNumberOfNodesPerElement = mLagrangeMesh->get_number_of_basis_per_element();

            // unity matrix
            Matrix< DDRMat > tEye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, 0.0 );
            for( uint i=0; i<tNumberOfNodesPerElement; ++i )
            {
                tEye( i, i ) = 1.0;
            }

            // calculate transposed Lagrange T-Matrix
            Matrix< DDRMat > tL( this->get_lagrange_matrix() );

            // loop over all elements
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to element
                auto tLagrangeElement = mLagrangeMesh->get_element( e );

                // FIXME : activate this flag
                //if ( tLagrangeElement->get_t_matrix_flag() )
                {
                    // get pointer to background element
                    auto tBackgroundElement = tLagrangeElement->get_background_element();

                    // initialize refinement Matrix
                    Matrix< DDRMat > tR( tEye );

                    while( ! tBackgroundElement->is_active( tBSplinePattern ) )
                    {
                        // right multiply refinement matrix
                        tR = this->get_refinement_matrix(
                                tBackgroundElement->get_child_index() ) * tR;

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

                    // transposed T-Matrix
                    Matrix< DDRMat > tT = tR * tL * tB;

                    // number of columns in T-Matrix
                    uint tNCols = tT.n_cols();

                    // epsilon to count T-Matrix
                    real tEpsilon = 1e-12;

                    // loop over all nodes of this element
                    for( uint k = 0; k<tNumberOfNodesPerElement; ++k  )
                    {
                        // pointer to node
                        auto tNode = tLagrangeElement->get_basis( k );

                        // test if node is flagged
                        if ( ! tNode->is_flagged() )
                        {
                            // initialize counter
                            uint tCount = 0;

                            // count number of nonzero entries
                            for( uint i=0; i<tNCols; ++i )
                            {
                                if ( std::abs( tT( k, i ) ) > tEpsilon )
                                {
                                    // increment counter
                                    ++tCount;
                                }
                            }

                            // reserve DOF cell
                            Cell< mtk::Vertex* > tNodeDOFs( tCount, nullptr );

                            // reserve matrix with coefficients
                            Matrix< DDRMat > tCoefficients( tCount, 1 );

                            // reset counter
                            tCount = 0;

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

                            // store the coefficients
                            tNode->set_weights( tCoefficients );

                            // store pointers to the DOFs
                            tNode->set_coefficients( tNodeDOFs );

                            // flag this node as processed
                            tNode->flag();
                        }
                    } // end loop over all nodes of this element
                } // end if element is flagged
            } // end loop over all elements
        }

//-------------------------------------------------------------------------------

/*        void
        T_Matrix::init_gauss_points()
        {
            // get number of points needed for exact mass matrix integration
            uint tNumberOfPoints =  ( mBSplineOrder > mLagrangeOrder ) ?
                    ( 4*mBSplineOrder - 1 ) : ( 4*mLagrangeOrder - 1 );

            // support points 1D
            Matrix< DDRMat > tXi( tNumberOfPoints, 1, 0 );

            // weights 1D
            Matrix< DDRMat > tW( tNumberOfPoints, 1, 0 );

            switch( tNumberOfPoints )
            {
                case( 1 ) :
                {
                    tXi( 0 ) = 0;
                    tW ( 0 ) = 1;
                    break;
                }
                case( 2 ):
                {
                    tXi( 0 ) = -std::sqrt( 1.0/3.0 );
                    tXi( 1 ) =  std::sqrt( 1.0/3.0 );
                    tW( 0 ) = 1.0;
                    tW( 1 ) = 1.0;
                    break;
                }
                default :
                {
                    // only need to evaluate half
                    uint tN = tNumberOfPoints / 2;

                    // Archimedes' constant
                    real tPi = 2*std::acos( 0 );

                    // epsilon for moris::real
                    long double tEpsilon = std::numeric_limits<real>::epsilon();


                    // loop over all points
                    for( uint k=1; k<=tN; ++k )
                    {
                        // guess xi
                        long double tXi1 = std::cos( tPi*(4*k-1)/(4*tNumberOfPoints+2) );

                        // initial value for old xi
                        long double tXi0= -2;

                        // Legendre polynomial
                        long double tP;

                        // derivative of Legendre polynomial
                        long double tdPdXi;

                        // perform Newton iteration
                        while( std::abs( tXi1-tXi0) > tEpsilon )
                        {
                            // evaluate Legendre polynomial
                            this->legendre( tNumberOfPoints, tXi1, tP, tdPdXi );

                            // shift xi
                            tXi0 = tXi1;

                            // perform Newton step
                            tXi1 -= tP/tdPdXi;
                        }

                        // save negative values
                        tXi( k-1 ) = - ( real ) tXi1;
                        tW( k-1 )  = ( real ) 2.0/( ( 1.0-std::pow( tXi1, 2 ) )*std::pow( tdPdXi, 2 ) );

                        // mirror positive values
                        tXi ( tNumberOfPoints-k ) = ( real ) tXi1;
                        tW(  tNumberOfPoints-k ) = tW( k-1 );

                    }

                    // special case for odd number of points
                    if ( tNumberOfPoints % 2 == 1 )
                    {
                        tXi( tN ) = 0.0;
                        tW ( tN ) = 0.0;
                        tW ( tN ) = 2.0 - sum( tW );
                    }
                }
            }

            // copy points to memory
            mGaussPoints = tXi;
            mGaussWeights = tW;
        }
//-------------------------------------------------------------------------------

        void
        T_Matrix::legendre(
                const uint         & aIndex,
                const long double  & aX,
                long double        & aP,
                long double        & adPdX ) const
        {
            switch ( aIndex )
            {
                case( 1 ) :
                {
                    aP    = 1.0;
                    adPdX = 0.0;
                    break;
                }
                case( 2 ) :
                {
                    aP    = aX;
                    adPdX = 1.0;
                    break;
                }
                default:
                {

                    // square of tXi
                    long double tX2 = std::pow( aX, 2 );

                    // values of polynomials
                    long double f1 = aX;               // P1
                    long double f0 = 0.5*( 3*tX2-1 );  // P2


                    if ( aIndex > 2 )
                    {

                        long double f2;
                        // calculate higher order polynomial
                        for( uint k=3; k<=aIndex; ++k )
                        {
                            // shift functions
                            f2 = f1;
                            f1 = f0;
                            f0 = ( ( ( 2*k-1 )*aX )*f1 - ( k-1 )*f2 )/
                                    ( ( long double ) k );
                        }
                    }

                    // copy output value
                    aP = f0;

                    // calculate derivative
                    adPdX =  aIndex* ( aX*f0-f1 ) / ( tX2 - 1 );

                    break;
                }
            }
        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::init_mass_matrices()
        {
            // ask settings for dimensions
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // ask meshes for order
            uint tNodesPerDirection  = mBSplineOrder + 1;
            uint tBasisPerDirection  = mLagrangeOrder + 1;

            // calculate number of Lagrange Nodes
            uint tNumberOfNodes = std::pow( tNodesPerDirection,
                                            tNumberOfDimensions );

            // calculate number of Control points
            uint tNumberOfBasis = std::pow(
                    tBasisPerDirection,
                    tNumberOfDimensions );

            // initialize matrices
            mBSplineMass.set_size( tNumberOfBasis, tNumberOfBasis, 0 );
            mLagrangeMass.set_size( tNumberOfNodes, tNumberOfBasis, 0 );

            // container for B-Spline shape function
            Matrix< DDRMat > tNB( tNumberOfBasis, 1 );

            // container for Lagrange shape
            Matrix< DDRMat > tNL( tNumberOfNodes, 1 );

            // get number of gauss points
            uint tNumberOfGaussPoints = mGaussWeights.length();

            if ( tNumberOfDimensions == 2 )
            {
                for( uint j=0; j< tNumberOfGaussPoints; ++j )
                {
                    for( uint i=0; i< tNumberOfGaussPoints; ++i )
                    {
                        // evaluate Lagrange functions
                        this->lagrange_shape(
                                mGaussPoints( i ),
                                mGaussPoints( j ),
                                tNL );

                        // evaluate B-Spline functions
                        this->b_spline_shape(
                                mGaussPoints( i ),
                                mGaussPoints( j ),
                                tNB );

                        // evaluate weight
                        real tW = mGaussWeights( i )*mGaussWeights( j );

                        mBSplineMass = mBSplineMass + tW * ( tNB * trans( tNB ) );
                        mLagrangeMass = mLagrangeMass + tW * ( tNB * trans( tNL ) );
                    }
                }
            }
            else if ( tNumberOfDimensions == 3 )
            {
                for( uint k=0; k< tNumberOfGaussPoints; ++k )
                {
                    for( uint j=0; j< tNumberOfGaussPoints; ++j )
                    {
                        for( uint i=0; i< tNumberOfGaussPoints; ++i )
                        {
                            // evaluate Lagrange functions
                            this->lagrange_shape(
                                    mGaussPoints( i ),
                                    mGaussPoints( j ),
                                    mGaussPoints( k ),
                                    tNL );

                            // evaluate B-Spline functions
                            this->b_spline_shape(
                                    mGaussPoints( i ),
                                    mGaussPoints( j ),
                                    mGaussPoints( k ),
                                    tNB );

                            // evaluate weight
                            real tW = mGaussWeights( i )
                                     *mGaussWeights( j )
                                     *mGaussWeights( k );

                            mBSplineMass  = mBSplineMass  + tW * ( tNB * trans( tNB ) );
                            mLagrangeMass = mLagrangeMass + tW * ( tNB * trans( tNL ) );
                        }
                    }
                }
            }
            //mBSplineMass.print("B");
            //mLagrangeMass.print("L");
            //Matrix< DDRMat > mBSplineMass;
            //Matrix< DDRMat > mLagrangeMass;
        } */

//-------------------------------------------------------------------------------

        void
        T_Matrix::init_lagrange_refinement_matrices()
        {

            // tidy up memory
            mLagrangeRefinementMatrix.clear();

            // get number of dimensions
            uint tNumberOfDimensions = mLagrangeParam.n_rows();

            // get number of nodes
            uint tNumberOfNodes      = mLagrangeParam.n_cols();

            // number of children
            uint tNumberOfChildren = std::pow( 2, tNumberOfDimensions );

            // initialize container
            Matrix< DDRMat > tEmpty( tNumberOfNodes, tNumberOfNodes, 0.0 );
            mLagrangeRefinementMatrix.resize( tNumberOfChildren,tEmpty );


            // matrix containing corner nodes
            Matrix< DDRMat > tCorners( tNumberOfChildren, tNumberOfDimensions );

            // shape function for "geometry"
            Matrix< DDRMat > tNGeo( 1, tNumberOfChildren );

            // shape function
             Matrix< DDRMat > tN( 1, tNumberOfNodes );

            // step 1: get parameter coordinates of child

            // matrix with parameter coordinates
            //Matrix< DDRMat > tXi( tNumberOfNodes, tNumberOfDimensions );

            // loop over all children
            for( uint c=0; c<tNumberOfChildren; ++c )
            {
                // get matrix with  corner nodes
                mGetCorners( c, tCorners );

                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    // evaluate shape function for "geometry"
                    mEvalNGeo( mLagrangeParam.cols( k, k ), tNGeo );

                    // get parameter coordinates
                    Matrix< DDRMat > tXi = tNGeo * tCorners;

                    // evaluate shape function
                    ( this->*mEvalN )(
                            tXi,
                            tN );

                    // copy result into matrix
                    mLagrangeRefinementMatrix( c ).rows( k, k ) = tN.rows( 0, 0 );

                }
                //mLagrangeRefinementMatrix( c ).print("T");
            }
        }

//------------------------------------------------------------------------------

        /**
         * returns the corner nodes of a child and dimension
         */
        void
        T_Matrix::get_child_corner_nodes_2d( const uint & aChildIndex, Matrix< DDRMat > & aXi )
        {
            switch ( aChildIndex )
            {
                case( 0 ) :
                {
                    aXi( 0, 0 ) = -1.0;
                    aXi( 1, 0 ) =  0.0;
                    aXi( 2, 0 ) =  0.0;
                    aXi( 3, 0 ) = -1.0;

                    aXi( 0, 1 ) = -1.0;
                    aXi( 1, 1 ) = -1.0;
                    aXi( 2, 1 ) =  0.0;
                    aXi( 3, 1 ) =  0.0;

                    break;
                }
                case( 1 ) :
                {
                    aXi( 0, 0 ) =  0.0;
                    aXi( 1, 0 ) =  1.0;
                    aXi( 2, 0 ) =  1.0;
                    aXi( 3, 0 ) =  0.0;

                    aXi( 0, 1 ) = -1.0;
                    aXi( 1, 1 ) = -1.0;
                    aXi( 2, 1 ) =  0.0;
                    aXi( 3, 1 ) =  0.0;

                    break;
                }

                case( 2 ) :
                {
                    aXi( 0, 0 ) = -1.0;
                    aXi( 1, 0 ) =  0.0;
                    aXi( 2, 0 ) =  0.0;
                    aXi( 3, 0 ) = -1.0;

                    aXi( 0, 1 ) =  0.0;
                    aXi( 1, 1 ) =  0.0;
                    aXi( 2, 1 ) =  1.0;
                    aXi( 3, 1 ) =  1.0;

                    break;
                }

                case( 3 ) :
                {
                    aXi( 0, 0 ) =  0.0;
                    aXi( 1, 0 ) =  1.0;
                    aXi( 2, 0 ) =  1.0;
                    aXi( 3, 0 ) =  0.0;

                    aXi( 0, 1 ) =  0.0;
                    aXi( 1, 1 ) =  0.0;
                    aXi( 2, 1 ) =  1.0;
                    aXi( 3, 1 ) =  1.0;

                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "invalid child index" );
                }
            }
        }

//-------------------------------------------------------------------------------

        /**
         * returns the corner nodes of a child and dimension
         */
        void
        T_Matrix::get_child_corner_nodes_3d( const uint & aChildIndex, Matrix< DDRMat > & aXi )
        {
            switch ( aChildIndex )
            {
                case( 0 ) :
                {
                    aXi( 0, 0 ) = -1.0;
                    aXi( 1, 0 ) =  0.0;
                    aXi( 2, 0 ) =  0.0;
                    aXi( 3, 0 ) = -1.0;
                    aXi( 4, 0 ) = -1.0;
                    aXi( 5, 0 ) =  0.0;
                    aXi( 6, 0 ) =  0.0;
                    aXi( 7, 0 ) = -1.0;

                    aXi( 0, 1 ) = -1.0;
                    aXi( 1, 1 ) = -1.0;
                    aXi( 2, 1 ) =  0.0;
                    aXi( 3, 1 ) =  0.0;
                    aXi( 4, 1 ) = -1.0;
                    aXi( 5, 1 ) = -1.0;
                    aXi( 6, 1 ) =  0.0;
                    aXi( 7, 1 ) =  0.0;

                    aXi( 0, 2 ) = -1.0;
                    aXi( 1, 2 ) = -1.0;
                    aXi( 2, 2 ) = -1.0;
                    aXi( 3, 2 ) = -1.0;
                    aXi( 4, 2 ) =  0.0;
                    aXi( 5, 2 ) =  0.0;
                    aXi( 6, 2 ) =  0.0;
                    aXi( 7, 2 ) =  0.0;

                    break;
                }
                case( 1 ) :
                {
                    aXi( 0, 0 ) =  0.0;
                    aXi( 1, 0 ) =  1.0;
                    aXi( 2, 0 ) =  1.0;
                    aXi( 3, 0 ) =  0.0;
                    aXi( 4, 0 ) =  0.0;
                    aXi( 5, 0 ) =  1.0;
                    aXi( 6, 0 ) =  1.0;
                    aXi( 7, 0 ) =  0.0;

                    aXi( 0, 1 ) = -1.0;
                    aXi( 1, 1 ) = -1.0;
                    aXi( 2, 1 ) =  0.0;
                    aXi( 3, 1 ) =  0.0;
                    aXi( 4, 1 ) = -1.0;
                    aXi( 5, 1 ) = -1.0;
                    aXi( 6, 1 ) =  0.0;
                    aXi( 7, 1 ) =  0.0;

                    aXi( 0, 2 ) = -1.0;
                    aXi( 1, 2 ) = -1.0;
                    aXi( 2, 2 ) = -1.0;
                    aXi( 3, 2 ) = -1.0;
                    aXi( 4, 2 ) =  0.0;
                    aXi( 5, 2 ) =  0.0;
                    aXi( 6, 2 ) =  0.0;
                    aXi( 7, 2 ) =  0.0;

                    break;
                }
                case( 2 ) :
                {
                    aXi( 0, 0 ) = -1.0;
                    aXi( 1, 0 ) =  0.0;
                    aXi( 2, 0 ) =  0.0;
                    aXi( 3, 0 ) = -1.0;
                    aXi( 4, 0 ) = -1.0;
                    aXi( 5, 0 ) =  0.0;
                    aXi( 6, 0 ) =  0.0;
                    aXi( 7, 0 ) = -1.0;

                    aXi( 0, 1 ) =  0.0;
                    aXi( 1, 1 ) =  0.0;
                    aXi( 2, 1 ) =  1.0;
                    aXi( 3, 1 ) =  1.0;
                    aXi( 4, 1 ) =  0.0;
                    aXi( 5, 1 ) =  0.0;
                    aXi( 6, 1 ) =  1.0;
                    aXi( 7, 1 ) =  1.0;

                    aXi( 0, 2 ) = -1.0;
                    aXi( 1, 2 ) = -1.0;
                    aXi( 2, 2 ) = -1.0;
                    aXi( 3, 2 ) = -1.0;
                    aXi( 4, 2 ) =  0.0;
                    aXi( 5, 2 ) =  0.0;
                    aXi( 6, 2 ) =  0.0;
                    aXi( 7, 2 ) =  0.0;

                    break;
                }
                case( 3 ) :
                {
                    aXi( 0, 0 ) =  0.0;
                    aXi( 1, 0 ) =  1.0;
                    aXi( 2, 0 ) =  1.0;
                    aXi( 3, 0 ) =  0.0;
                    aXi( 4, 0 ) =  0.0;
                    aXi( 5, 0 ) =  1.0;
                    aXi( 6, 0 ) =  1.0;
                    aXi( 7, 0 ) =  0.0;

                    aXi( 0, 1 ) =  0.0;
                    aXi( 1, 1 ) =  0.0;
                    aXi( 2, 1 ) =  1.0;
                    aXi( 3, 1 ) =  1.0;
                    aXi( 4, 1 ) =  0.0;
                    aXi( 5, 1 ) =  0.0;
                    aXi( 6, 1 ) =  1.0;
                    aXi( 7, 1 ) =  1.0;

                    aXi( 0, 2 ) = -1.0;
                    aXi( 1, 2 ) = -1.0;
                    aXi( 2, 2 ) = -1.0;
                    aXi( 3, 2 ) = -1.0;
                    aXi( 4, 2 ) =  0.0;
                    aXi( 5, 2 ) =  0.0;
                    aXi( 6, 2 ) =  0.0;
                    aXi( 7, 2 ) =  0.0;

                    break;
                }
                case( 4 ) :
                {
                    aXi( 0, 0 ) = -1.0;
                    aXi( 1, 0 ) =  0.0;
                    aXi( 2, 0 ) =  0.0;
                    aXi( 3, 0 ) = -1.0;
                    aXi( 4, 0 ) = -1.0;
                    aXi( 5, 0 ) =  0.0;
                    aXi( 6, 0 ) =  0.0;
                    aXi( 7, 0 ) = -1.0;

                    aXi( 0, 1 ) = -1.0;
                    aXi( 1, 1 ) = -1.0;
                    aXi( 2, 1 ) =  0.0;
                    aXi( 3, 1 ) =  0.0;
                    aXi( 4, 1 ) = -1.0;
                    aXi( 5, 1 ) = -1.0;
                    aXi( 6, 1 ) =  0.0;
                    aXi( 7, 1 ) =  0.0;

                    aXi( 0, 2 ) =  0.0;
                    aXi( 1, 2 ) =  0.0;
                    aXi( 2, 2 ) =  0.0;
                    aXi( 3, 2 ) =  0.0;
                    aXi( 4, 2 ) =  1.0;
                    aXi( 5, 2 ) =  1.0;
                    aXi( 6, 2 ) =  1.0;
                    aXi( 7, 2 ) =  1.0;

                    break;
                }
                case( 5 ) :
                {
                    aXi( 0, 0 ) =  0.0;
                    aXi( 1, 0 ) =  1.0;
                    aXi( 2, 0 ) =  1.0;
                    aXi( 3, 0 ) =  0.0;
                    aXi( 4, 0 ) =  0.0;
                    aXi( 5, 0 ) =  1.0;
                    aXi( 6, 0 ) =  1.0;
                    aXi( 7, 0 ) =  0.0;

                    aXi( 0, 1 ) = -1.0;
                    aXi( 1, 1 ) = -1.0;
                    aXi( 2, 1 ) =  0.0;
                    aXi( 3, 1 ) =  0.0;
                    aXi( 4, 1 ) = -1.0;
                    aXi( 5, 1 ) = -1.0;
                    aXi( 6, 1 ) =  0.0;
                    aXi( 7, 1 ) =  0.0;

                    aXi( 0, 2 ) =  0.0;
                    aXi( 1, 2 ) =  0.0;
                    aXi( 2, 2 ) =  0.0;
                    aXi( 3, 2 ) =  0.0;
                    aXi( 4, 2 ) =  1.0;
                    aXi( 5, 2 ) =  1.0;
                    aXi( 6, 2 ) =  1.0;
                    aXi( 7, 2 ) =  1.0;

                    break;
                }
                case( 6 ) :
                {
                    aXi( 0, 0 ) = -1.0;
                    aXi( 1, 0 ) =  0.0;
                    aXi( 2, 0 ) =  0.0;
                    aXi( 3, 0 ) = -1.0;
                    aXi( 4, 0 ) = -1.0;
                    aXi( 5, 0 ) =  0.0;
                    aXi( 6, 0 ) =  0.0;
                    aXi( 7, 0 ) = -1.0;

                    aXi( 0, 1 ) =  0.0;
                    aXi( 1, 1 ) =  0.0;
                    aXi( 2, 1 ) =  1.0;
                    aXi( 3, 1 ) =  1.0;
                    aXi( 4, 1 ) =  0.0;
                    aXi( 5, 1 ) =  0.0;
                    aXi( 6, 1 ) =  1.0;
                    aXi( 7, 1 ) =  1.0;

                    aXi( 0, 2 ) =  0.0;
                    aXi( 1, 2 ) =  0.0;
                    aXi( 2, 2 ) =  0.0;
                    aXi( 3, 2 ) =  0.0;
                    aXi( 4, 2 ) =  1.0;
                    aXi( 5, 2 ) =  1.0;
                    aXi( 6, 2 ) =  1.0;
                    aXi( 7, 2 ) =  1.0;

                    break;
                }
                case( 7 ) :
                {
                    aXi( 0, 0 ) =  0.0;
                    aXi( 1, 0 ) =  1.0;
                    aXi( 2, 0 ) =  1.0;
                    aXi( 3, 0 ) =  0.0;
                    aXi( 4, 0 ) =  0.0;
                    aXi( 5, 0 ) =  1.0;
                    aXi( 6, 0 ) =  1.0;
                    aXi( 7, 0 ) =  0.0;

                    aXi( 0, 1 ) =  0.0;
                    aXi( 1, 1 ) =  0.0;
                    aXi( 2, 1 ) =  1.0;
                    aXi( 3, 1 ) =  1.0;
                    aXi( 4, 1 ) =  0.0;
                    aXi( 5, 1 ) =  0.0;
                    aXi( 6, 1 ) =  1.0;
                    aXi( 7, 1 ) =  1.0;

                    aXi( 0, 2 ) =  0.0;
                    aXi( 1, 2 ) =  0.0;
                    aXi( 2, 2 ) =  0.0;
                    aXi( 3, 2 ) =  0.0;
                    aXi( 4, 2 ) =  1.0;
                    aXi( 5, 2 ) =  1.0;
                    aXi( 6, 2 ) =  1.0;
                    aXi( 7, 2 ) =  1.0;

                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "invalid child index" );
                }
            }
        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::N_quad4( const Matrix< DDRMat > & aXi, Matrix< DDRMat > & aN )
        {
            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // populate matrix with values
            aN( 0 ) = ( ( 1.0 - xi ) * ( 1.0 - eta ) ) * 0.25;
            aN( 1 ) = ( ( 1.0 + xi ) * ( 1.0 - eta ) ) * 0.25;
            aN( 2 ) = ( ( 1.0 + xi ) * ( 1.0 + eta ) ) * 0.25;
            aN( 3 ) = ( ( 1.0 - xi ) * ( 1.0 + eta ) ) * 0.25;
        }

//-------------------------------------------------------------------------------

        void
        T_Matrix::N_hex8( const Matrix< DDRMat > & aXi, Matrix< DDRMat > & aN )
        {
            // unpack xi and eta from input vector
            auto    xi = aXi( 0 );
            auto   eta = aXi( 1 );
            auto  zeta = aXi( 2 );

            // populate output matrix
            aN( 0 ) =  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 1 ) =    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 2 ) =  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 3 ) =    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 4 ) =    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 5 ) =  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 6 ) =    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 7 ) =  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
        }

//-------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

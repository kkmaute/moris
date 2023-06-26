/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_T_MATRIX_HPP_
#define SRC_HMR_CL_HMR_T_MATRIX_HPP_

#include "cl_HMR_T_Matrix_Base.hpp"
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "cl_Cell.hpp" //CNT/src
#include "fn_eye.hpp"

namespace moris::hmr
{
    /**
     * T-matrix generator class
     *
     * @tparam N Number of dimensions
     */
    template< uint N >
    class T_Matrix : public T_Matrix_Base
    {
    public:

        /**
         * Constructor
         *
         * @param aLagrangeMesh Lagrange mesh pointer
         * @param aBSplineMesh B-spline Mesh pointer
         * @param aTruncate Whether or not to truncate B-splines
         */
        T_Matrix(
                Lagrange_Mesh_Base* aLagrangeMesh,
                BSpline_Mesh_Base*  aBSplineMesh = nullptr,
                bool                aTruncate = true )
                : T_Matrix_Base( aLagrangeMesh, aBSplineMesh, aTruncate )
        {
            // Initializations
            this->init_lagrange_parameter_coordinates();
            this->init_lagrange_refinement_matrices();
            this->init_lagrange_change_order_matrices();

            // If B-spline mesh is given
            if ( aBSplineMesh != nullptr )
            {
                this->init_basis_index();
                this->init_unity_matrix();
                this->init_child_matrices();
                this->init_truncation_weights();
                this->init_lagrange_matrix();
            }
        }

//-------------------------------------------------------------------------------

        // destructor
        ~T_Matrix()
        {
        }

        /**
         * @brief evaluate the L2 projection weights converting extended basis to root basis
         *
         * @param aRootBsplineElement
         * @param aExtendedBsplineElement
         * @param aRootBsplineBasis
         * @param aExtendedBsplineBasis
         * @param aWeights
         */
        void evaluate_L2_projection(
                Element*                                    aRootBsplineElement,
                Element*                                    aExtendedBsplineElement,
                moris::Cell< moris::Cell< mtk::Vertex* > >& aRootBsplineBasis,
                moris::Cell< mtk::Vertex* >&                aExtendedBsplineBasis,
                moris::Cell< Matrix< DDRMat > >&            aWeights )
        {
            // Background_Element_Base* aRootBSpBackgroundElement = aRootBsplineElement->get_background_element();
            // Background_Element_Base* aExtendedBSpBackgroundElement = aExtendedBsplineElement->get_background_element();

            const luint* tRootIJK     = aRootBsplineElement->get_ijk();
            const luint* tExtendedIJK = aExtendedBsplineElement->get_ijk();

            // number of bspline coefficients per direction
            uint tNodesPerDirection = mBSplineOrder + 1;

            // initialize 1D matrices to find the projection matrix in 1D
            moris::Cell< Matrix< DDRMat > > tMatrices1D( N, Matrix< DDRMat >( tNodesPerDirection, tNodesPerDirection, MORIS_REAL_MAX ) );

            // loop over the dimensions and create the i,j,k 1D matrices
            for ( uint iDim = 0; iDim < N; iDim++ )
            {
                real tShift = tRootIJK[ iDim ] < tExtendedIJK[ iDim ] ? real( -tRootIJK[ iDim ] + tExtendedIJK[ iDim ] ) : -real( -tExtendedIJK[ iDim ] + tRootIJK[ iDim ] );

                // find the shift in each direction
                this->get_extension_matrix_1d( tShift, tMatrices1D( iDim ) );
            }

            // calculate number of basis per element
            uint tNumberOfBasis = std::pow( tNodesPerDirection, N );

            // initialize the projection matrix with the correct size
            Matrix< DDRMat > tL2ProjectionMatrix( tNumberOfBasis, tNumberOfBasis );

            // Apply a tensor product to get the final weights
            if ( N == 2 )
            {
                uint b = 0;
                for ( uint l = 0; l < tNodesPerDirection; ++l )
                {
                    for ( uint k = 0; k < tNodesPerDirection; ++k )
                    {
                        uint a = 0;
                        for ( uint j = 0; j < tNodesPerDirection; ++j )
                        {
                            for ( uint i = 0; i < tNodesPerDirection; ++i )
                            {
                                tL2ProjectionMatrix( mBasisIndex( a ), mBasisIndex( b ) ) = tMatrices1D( 0 )( i, k ) * tMatrices1D( 1 )( j, l );
                                ++a;
                            }
                        }
                        ++b;
                    }
                }
            }
            else if ( N == 3 )
            {
                uint b = 0;
                for ( uint p = 0; p < tNodesPerDirection; ++p )
                {
                    for ( uint q = 0; q < tNodesPerDirection; ++q )
                    {
                        for ( uint l = 0; l < tNodesPerDirection; ++l )
                        {
                            uint a = 0;
                            for ( uint k = 0; k < tNodesPerDirection; ++k )
                            {
                                for ( uint j = 0; j < tNodesPerDirection; ++j )
                                {
                                    for ( uint i = 0; i < tNodesPerDirection; ++i )
                                    {
                                        tL2ProjectionMatrix( mBasisIndex( a ), mBasisIndex( b ) ) = tMatrices1D( 0 )( i, l ) * tMatrices1D( 1 )( j, q ) * tMatrices1D( 2 )( k, p );
                                        ++a;
                                    }
                                }
                            }
                            ++b;
                        }
                    }
                }
            }

            // initialize the basis for the root cell and extended cell
            moris::Cell< Basis* > tRootBasis;
            moris::Cell< Basis* > tExtendedBasis;

            //reserve enough memory for each of them
            tRootBasis.reserve(tNumberOfBasis);
            tExtendedBasis.reserve(tNumberOfBasis);

            // fill out the cell data
            for ( uint i = 0; i < tNumberOfBasis; i++ )
            {
                // get the basis
                Basis* tBasis = aRootBsplineElement->get_basis( i );

                //if it is active add it to the cell
                if ( tBasis->is_active() )
                {
                    tRootBasis.push_back( tBasis );
                }

                // get the basis
                tBasis = aExtendedBsplineElement->get_basis( i );

                //if it is active add it to the cell
                if ( tBasis->is_active() )
                {
                    tExtendedBasis.push_back( tBasis);
                }
            }

            // resize the output cells
            aWeights.resize(tNumberOfBasis );
            aRootBsplineBasis.resize(tNumberOfBasis );
            aExtendedBsplineBasis.resize(tNumberOfBasis );

            // loop over the basis and eliminate the zero values
            for ( uint iExtendedBasisIndex = 0; iExtendedBasisIndex < tNumberOfBasis; iExtendedBasisIndex++ )
            {
                // set size for each of the extended basis
                aExtendedBsplineBasis(iExtendedBasisIndex ) = tExtendedBasis(iExtendedBasisIndex );
                aWeights(iExtendedBasisIndex ).set_size(tNumberOfBasis, 1 );
                aRootBsplineBasis(iExtendedBasisIndex ).resize(tNumberOfBasis );

                // find all the basis and weights of the root that hve non-zero values
                uint tNonzeroCount = 0;
                for ( uint iRootBasisIndex = 0; iRootBasisIndex < tNumberOfBasis; iRootBasisIndex++ )
                {
                    // if the wieght is non-zero add it to the cell
                    if ( std::abs( tL2ProjectionMatrix( iExtendedBasisIndex, iRootBasisIndex ) ) > MORIS_REAL_EPS )
                    {
                        // assign the output values
                        aWeights(iExtendedBasisIndex )(tNonzeroCount )          = tL2ProjectionMatrix(iExtendedBasisIndex, iRootBasisIndex );
                        aRootBsplineBasis(iExtendedBasisIndex )(tNonzeroCount ) = tRootBasis(iRootBasisIndex );

                        //increment the non-zero count
                        tNonzeroCount++;
                    }
                }

                // resize the output cells to the correct non-zero size
                aWeights(iExtendedBasisIndex ).resize(tNonzeroCount, 1 );
                aRootBsplineBasis(iExtendedBasisIndex ).resize(tNonzeroCount );
            }
        }

    private:

        /**
         * determines the sorting order of the nodes
         */
        void init_basis_index()
        {
            // Create background element
            Background_Element_Base* tBackElement = this->create_background_element();

            // create a prototype of a B-Spline Element
            Element* tElement = mBSplineMesh->create_element( tBackElement );

            // B-Spline order
            mBSplineOrder = mBSplineMesh->get_order();

            // number of bspline coefficients per direction
            uint tNodesPerDirection = mBSplineOrder + 1;

            // calculate number of basis per element
            uint tNumberOfBasis = static_cast<uint>( std::pow( tNodesPerDirection, N ) );

            // initialize index matrix
            mBasisIndex.set_size( tNumberOfBasis, 1 );
            mBSplineIJK.set_size( N, tNumberOfBasis );

            // loop over all basis
            if ( N == 2 )
            {
                // container for ijk position of basis
                luint tIJ[ 2 ];
                for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBasis; iBasisIndex++ )
                {
                    // get position from element
                    tElement->get_ijk_of_basis( iBasisIndex, tIJ );

                    // calculate index in matrix
                    uint tIndex = tIJ[ 0 ] + tIJ[ 1 ] * tNodesPerDirection;

                    mBasisIndex( tIndex ) = iBasisIndex;
                    mBSplineIJK( 0, iBasisIndex )   = tIJ[ 0 ];
                    mBSplineIJK( 1, iBasisIndex )   = tIJ[ 1 ];
                }
            }
            else if ( N == 3 )
            {
                // container for ijk position of basis
                luint tIJK[ 3 ];
                for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBasis; iBasisIndex++ )
                {
                    // get position from element
                    tElement->get_ijk_of_basis( iBasisIndex, tIJK );

                    // calculate index in matrix
                    uint tIndex =
                            tIJK[ 0 ] + tNodesPerDirection * ( tIJK[ 1 ] + tIJK[ 2 ] * tNodesPerDirection );

                    mBasisIndex( tIndex ) = iBasisIndex;

                    mBSplineIJK( 0, iBasisIndex ) = tIJK[ 0 ];
                    mBSplineIJK( 1, iBasisIndex ) = tIJK[ 1 ];
                    mBSplineIJK( 2, iBasisIndex ) = tIJK[ 2 ];
                }
            }

            // tidy up
            delete tElement;
            delete tBackElement;
        }

//-------------------------------------------------------------------------------

        /**
         * initializes the container for the unity matrix
         */
        void init_unity_matrix()
        {
            // get number of basis per element
            uint tNumberOfBasis = static_cast<uint>(
                    std::pow( mBSplineMesh->get_order() + 1, N ) );

            // Set unity and zero matrices
            eye( tNumberOfBasis, tNumberOfBasis, mEye );
            mZero.set_size( tNumberOfBasis, 1, 0.0 );
        }

//-------------------------------------------------------------------------------
        /**
         * pre-calculates the child relation matrices
         */
        void init_child_matrices()
        {
            // get order of mesh
            uint tOrder = mBSplineMesh->get_order();

            // create temporary matrix
            Matrix< DDRMat > tFactors( tOrder + 1, tOrder + 2, 0.0 );

            // number of bspline coefficients per direction
            uint tNumCoefficients = tOrder + 1;

            // weight factor
            real tWeight = 1.0 / std::pow( 2, tOrder );

            for ( uint iCoefficient = 0; iCoefficient <= tNumCoefficients; iCoefficient++ )
            {
                for ( uint iOrder = 0; iOrder <= tOrder; iOrder++ )
                {
                    uint k = tOrder - 2 * iOrder + iCoefficient;
                    if ( k <= tNumCoefficients )
                    {
                        tFactors( iOrder, iCoefficient ) = tWeight * nchoosek( tNumCoefficients, k );
                    }
                }
            }

            // left and right matrices
            Matrix< DDRMat > TL( tNumCoefficients, tNumCoefficients, 0.0 );
            Matrix< DDRMat > TR( tNumCoefficients, tNumCoefficients, 0.0 );

            // Fill matrices
            for ( uint iOrder = 0; iOrder <= tOrder; iOrder++ )
            {
                TL.set_column( iOrder, tFactors.get_column( iOrder ) );
                TR.set_column( iOrder, tFactors.get_column( iOrder + 1 ) );
            }

            // determine number of children
            uint tNumberOfChildren = static_cast<uint>( std::pow( 2, N ) );

            // determine number of basis per element
            uint tNumberOfBasis = static_cast<uint>( std::pow( tOrder + 1, N ) );

            // empty matrix
            Matrix< DDRMat > tEmpty( tNumberOfBasis, tNumberOfBasis, 0.0 );

            // container for child relation matrices ( transposed! )
            mChild.resize( tNumberOfChildren, tEmpty );

            // populate child matrices.
            this->populate_child_matrices( TL, TR );

            // Child multiplication
            this->child_multiplication();
        }

//-------------------------------------------------------------------------------

        void child_multiplication()
        {
            // determine number of children
            uint tNumberOfChildren = static_cast<uint>( std::pow( 2, N ) );
            mChildMultiplied.resize( tNumberOfChildren );

            for ( uint Ik = 0; Ik < tNumberOfChildren; ++Ik )
            {
                mChildMultiplied( Ik ).resize( tNumberOfChildren );
            }

            // multiply child matrices
            for ( uint iChildRow = 0; iChildRow < tNumberOfChildren; iChildRow++ )
            {
                for ( uint iChildCol = 0; iChildCol < tNumberOfChildren; iChildCol++ )
                {
                    mChildMultiplied( iChildRow )( iChildCol ) = mChild( iChildRow ) * mChild( iChildCol );
                }
            }
        }

//------------------------------------------------------------------------------

        void init_truncation_weights()
        {
            // get order of mesh
            uint tOrder = mBSplineMesh->get_order();

            // number of children per direction
            uint tNumberOfChildren = tOrder + 2;

            // matrix containing 1D weights
            Matrix< DDRMat > tWeights( tNumberOfChildren, 1 );

            // scale factor for 1D weights
            real tScale = 1.0 / ( (real)std::pow( 2, tOrder ) );

            // calculate 1D weights
            for ( uint iChildIndex = 0; iChildIndex < tNumberOfChildren; iChildIndex++ )
            {
                tWeights( iChildIndex ) = tScale * nchoosek( tOrder + 1, iChildIndex );
            }

            // allocate weights
            mTruncationWeights.set_size( static_cast<uint>( std::pow( tNumberOfChildren, N ) ), 1 );

            // Evaluate weights
            evaluate_truncation_weights( tWeights );
        }

//------------------------------------------------------------------------------

        void init_lagrange_parameter_coordinates()
        {
            // create background element for reference
            Background_Element_Base* tBackElement = this->create_background_element();

            // number of nodes per direction
            uint tNodesPerDirection = mLagrangeMesh->get_order() + 1;

            // number of nodes
            mNumberOfNodes = static_cast<uint>( std::pow( tNodesPerDirection, N ) );

            // create a Lagrange element
            Element* tElement = mLagrangeMesh->create_element( tBackElement );

            // assign memory for parameter coordinates
            mLagrangeParam.set_size( N, mNumberOfNodes );

            // scaling factor
            real tScale = 1.0 / ( (real)mLagrangeMesh->get_order() );

            // container for ijk position of basis
            luint tIJK[ 3 ];

            // ijk positions for reference Lagrange element
            mLagrangeIJK.set_size( N, mNumberOfNodes );

            for ( uint iNodeIndex = 0; iNodeIndex < mNumberOfNodes; iNodeIndex++ )
            {
                // get position from element
                tElement->get_ijk_of_basis( iNodeIndex, tIJK );

                // save coordinate into memory
                for ( uint iDimensionIndex = 0; iDimensionIndex < N; iDimensionIndex++ )
                {
                    // fill in node ijk positions in element
                    mLagrangeParam( iDimensionIndex, iNodeIndex ) = 2 * tScale * tIJK[ iDimensionIndex ] - 1.0;

                    // fill in nodal natural coordinates for this element
                    mLagrangeIJK( iDimensionIndex, iNodeIndex ) = tIJK[ iDimensionIndex ];
                }
            }

            // tidy up
            delete tElement;
            delete tBackElement;
        }

        //------------------------------------------------------------------------------

        void init_modified_lagrange_parameter_coordinates(
                Element* aLagrangeElement,
                Element* aBSplineElement ) override
        {
            // number of nodes per direction
            uint tNodesPerDirection = mLagrangeMesh->get_order() + 1;

            // number of nodes
            mNumberOfNodes = std::pow( tNodesPerDirection, N );

            // create background element for reference
            Background_Element_Base* tBackElement = this->create_background_element();

            // create a Lagrange element
            Element* tReferenceElement = mLagrangeMesh->create_element( tBackElement );

            // assign memory for parameter coordinates
            mLagrangeParamModified.set_size( mLagrangeParam.n_rows(), mLagrangeParam.n_cols() );

            // scaling factor
            real tScale = 1.0 / ( (real)mLagrangeMesh->get_order() );

            // container for ijk position of basis
            luint tIJK[ 3 ];
            const luint* tIJKBSpline  = aBSplineElement->get_ijk();
            const luint* tIJKLagrange = aLagrangeElement->get_ijk();

            for ( uint iNodeIndex = 0; iNodeIndex < mNumberOfNodes ; iNodeIndex++ )
            {
                // get position from element
                tReferenceElement->get_ijk_of_basis(iNodeIndex, tIJK );

                // save coordinate into memory
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    // Get modified ijk value based on difference between Lagrange and B-spline element ijk
                    moris_index tIJKValue = tIJK[ iDimension ] + tIJKLagrange[iDimension] - tIJKBSpline[iDimension];

                    // fill in node ijk positions in element
                    mLagrangeParamModified(iDimension, iNodeIndex ) = 2 * tScale * tIJKValue - 1.0;
                }
            }
        }

        //------------------------------------------------------------------------------

        /**
         * calculates the matrix that converts B-Spline DOFs per element to Lagrange DOFs.
         */
        void init_lagrange_matrix()
        {
            // get number of basis per element of B-Spline mesh
            uint tNumberOfBasis = mBSplineIJK.n_cols();

            // get number of Lagrange nodes
            uint tNumberOfNodes = mLagrangeParam.n_cols();

            // initialize T-Matrix for B-Spline to Lagrange conversion
            mTMatrixLagrange.set_size( tNumberOfNodes, tNumberOfBasis, 1 );

            // get order of B-Spline mesh
            uint tOrder = mBSplineMesh->get_order();

            // loop over all Lagrange nodes
            for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; iNodeIndex++ )
            {
                // loop over all B-Spline Basis
                for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBasis; iBasisIndex++ )
                {
                    // loop over all dimensions
                    for ( uint iDimensionIndex = 0; iDimensionIndex < N; iDimensionIndex++ )
                    {
                        mTMatrixLagrange( iNodeIndex, iBasisIndex ) *= T_Matrix::b_spline_shape_1d(
                                tOrder,
                                mBSplineIJK( iDimensionIndex, iBasisIndex ),
                                mLagrangeParam( iDimensionIndex, iNodeIndex ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void init_lagrange_refinement_matrices()
        {
            // tidy up memory
            mLagrangeRefinementMatrix.clear();

            // get number of nodes
            uint tNumberOfNodes = mLagrangeParam.n_cols();

            // number of children
            uint tNumberOfChildren = static_cast<uint>( std::pow( 2, N ) );

            // initialize container
            Matrix< DDRMat > tEmpty( tNumberOfNodes, tNumberOfNodes, 0.0 );
            mLagrangeRefinementMatrix.resize( tNumberOfChildren, tEmpty );

            // matrix containing corner nodes
            Matrix< DDRMat > tCorners( tNumberOfChildren, N );

            // shape function for "geometry"
            Matrix< DDRMat > tNGeo( 1, tNumberOfChildren );

            // shape function
            Matrix< DDRMat > tN( 1, tNumberOfNodes );

            // step 1: get parameter coordinates of child

            // matrix with parameter coordinates
            // Matrix< DDRMat > tXi( tNumberOfNodes, N );

            // loop over all children
            for ( uint iChildIndex = 0; iChildIndex < tNumberOfChildren; iChildIndex++ )
            {
                // get matrix with  corner nodes
                get_child_corner_nodes( iChildIndex, tCorners );

                for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; iNodeIndex++ )
                {
                    // evaluate shape function for "geometry"
                    evaluate_geometry_interpolation( mLagrangeParam.get_column( iNodeIndex ), tNGeo );

                    // get parameter coordinates
                    Matrix< DDRMat > tXi = tNGeo * tCorners;

                    // evaluate shape function
                    evaluate_shape_function( tXi, tN );

                    // copy result into matrix
                    mLagrangeRefinementMatrix( iChildIndex ).set_row( iNodeIndex, tN.get_row( 0 ) );
                }
            }
        }

//-------------------------------------------------------------------------------

        void init_lagrange_change_order_matrices()
        {
            // empty matrix
            Matrix< DDRMat > tEmpty;

            // tidy up memory
            mLagrangeChangeOrderMatrix.clear();
            mLagrangeChangeOrderMatrix.push_back( tEmpty );

            // get number of nodes
            uint tNumberOfNodesThisMesh = mLagrangeParam.n_cols();

            for ( uint tOrder = 1; tOrder <= 3; ++tOrder )
            {
                // grab the parameter coordinates
                Matrix< DDRMat > tXi = get_supporting_points( tOrder );

                // get the number of nodes of the other mesh
                uint tNumberOfNodesOtherMesh = tXi.n_cols();

                // pointer to T-Matrix
                Matrix< DDRMat > tT(tNumberOfNodesOtherMesh, tNumberOfNodesThisMesh );

                // the shape function
                Matrix< DDRMat > tN(1, tNumberOfNodesThisMesh );

                // Point coordinate matrix
                Matrix< DDRMat > tPoint( N, 1 );

                // populate matrix
                for ( uint iNodeOtherMesh = 0; iNodeOtherMesh < tNumberOfNodesOtherMesh; iNodeOtherMesh++ )
                {
                    // copy coordinates into point TODO can we be more efficient?
                    for ( uint iDimension = 0; iDimension < N; iDimension++ )
                    {
                        tPoint( iDimension ) = tXi(iDimension, iNodeOtherMesh );
                    }

                    // evaluate shape function
                    evaluate_shape_function( tPoint,tN );

                    for ( uint iNodeThisMesh = 0; iNodeThisMesh < tNumberOfNodesThisMesh; iNodeThisMesh++ )
                    {
                        tT(iNodeOtherMesh, iNodeThisMesh ) = tN(iNodeThisMesh );
                    }
                }
                mLagrangeChangeOrderMatrix.push_back( tT );
            }
        }

//-------------------------------------------------------------------------------

        Background_Element_Base* create_background_element()
        {
            // create a prototype for a background element
            luint tIJK[ N ] = { 0 };
            return new Background_Element< N >(
                    nullptr,
                    0,
                    tIJK,
                    0,
                    0,
                    0,
                    gNoProcOwner );
        }

//------------------------------------------------------------------------------

        /**
         * 1D shape function
         */
        static real b_spline_shape_1d(
                uint aOrder,
                uint aK,
                real aXi )
        {
            // max number of entries in lookup table
            uint tSteps = 2 * ( aOrder + 1 );

            // temporary matrix that contains B-Spline segments
            Matrix< DDRMat > tDeltaXi( tSteps, 1, 0 );
            for ( uint i = 0; i < tSteps; ++i )
            {
                tDeltaXi( i ) = ( ( (real)i ) - ( (real)aOrder ) ) * 2.0 - 1.0;
            }

            // temporary matrix that contains evaluated values
            Matrix< DDRMat > tN( aOrder + 1, 1, 0 );

            // initialize zero order values
            for ( uint i = 0; i <= aOrder; ++i )
            {
                if ( tDeltaXi( i + aK ) <= aXi && aXi < tDeltaXi( i + aK + 1 ) )
                {
                    tN( i ) = 1.0;
                }
            }

            // loop over all orders
            for ( uint r = 1; r <= aOrder; ++r )
            {
                // copy values of tN into old matrix
                Matrix< DDRMat > tNold( tN );

                // loop over all contributions
                for ( uint i = 0; i <= aOrder - r; ++i )
                {
                    // help values
                    real tA = aXi - tDeltaXi( i + aK );
                    real tB = tDeltaXi( i + aK + r + 1 ) - aXi;

                    tN( i ) = 0.5 * ( tA * tNold( i ) + tB * ( tNold( i + 1 ) ) ) / ( (real)r );
                }
            }

            // first value in entry is shape value
            return tN( 0 );
        }

        //------------------------------------------------------------------------------

        static real b_spline_shape_1d_extended(
                uint aOrder,
                uint aK,
                real aXi )
        {
            switch ( aOrder )
            {
                // linear interpolation
                case 1:
                {
                    // local ordering of basis function
                    switch ( aK )
                    {
                        case 0:
                        {
                            return 0.5 * ( 1.0 - aXi );
                        }
                        case 1:
                        {
                            return 0.5 * ( 1.0 + aXi );
                        }
                        default:
                        {
                            MORIS_ERROR( false, "The specified local basis %u is not implemented", aK );
                            return 0.0;
                        }
                    }
                }

                default:
                {
                    MORIS_ERROR( false, "The specified order %u is not implemented", aOrder );
                    return 0.0;
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * 1D Lagrange Shape function
         *
         * @param aBasisNumber Basis number
         * @param aXi
         * @return Evaluated shape function
         */
        real lagrange_shape_1d(
                uint aBasisNumber,
                real aXi ) const
        {
            // use horner scheme to evaluate 1D Lagrange function
            real aResult = 0.0;
            for ( uint i = 0; i < mLagrangeOrder; i++ )
            {
                aResult = ( aResult + mLagrangeCoefficients( i, aBasisNumber ) ) * aXi;
            }

            return aResult + mLagrangeCoefficients( mLagrangeOrder, aBasisNumber );
        }

        //------------------------------------------------------------------------------

        void recompute_lagrange_matrix() override
        {
            // get number of basis per element of B-Spline mesh
            uint tNumberOfBasis = mBSplineIJK.n_cols();

            // get number of Lagrange nodes
            uint tNumberOfNodes = mLagrangeParam.n_cols();

            // initialize T-Matrix for B-Spline to Lagrange conversion
            mTMatrixLagrangeModified.set_size( tNumberOfNodes, tNumberOfBasis, 1 );

            // get order of B-Spline mesh
            uint tOrder = mBSplineMesh->get_order();

            // loop over all Lagrange nodes
            for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; iNodeIndex++ )
            {
                // loop over all B-Spline Basis
                for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBasis; iBasisIndex++ )
                {
                    // loop over all dimensions
                    for ( uint iDimension = 0; iDimension < N; iDimension++ )
                    {
                        mTMatrixLagrangeModified( iNodeIndex, iBasisIndex ) *= this->b_spline_shape_1d_extended(
                                tOrder,
                                mBSplineIJK( iDimension, iBasisIndex ),
                                mLagrangeParamModified( iDimension, iNodeIndex ) );
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------

        void get_extension_matrix_1d( real const & aShift, Matrix< DDRMat >& aExtensionMatrix )
        {
            switch ( mBSplineOrder )
            {
                case 1:
                    aExtensionMatrix = { { 1.0 - aShift, aShift }, { -aShift, 1.0 + aShift } };
                    break;
                case 2:
                    aExtensionMatrix = { { 0.5 * ( aShift - 2.0 ) * ( aShift - 1.0 ), aShift * ( -aShift + 2.0 ), 0.5 * aShift * ( aShift - 1.0 ) },    //
                            { 0.5 * aShift * ( aShift - 1.0 ), -( aShift - 1.0 ) * ( aShift + 1.0 ), 0.5 * aShift * ( aShift + 1 ) },                       //
                            { 0.5 * aShift * ( aShift + 1.0 ), -aShift * ( aShift + 2.0 ), 0.5 * ( aShift + 1 ) * ( aShift + 2 ) } };
                    break;

                default:
                    MORIS_ERROR( false, "B-spline order not known for extension matrix" );
                    break;
            }
        }

        //-------------------------------------------------------------------------------

        /**
         * Evaluates a (Lagrange) geometry interpolation at a given point
         *
         * @param aXi Local coordinates
         * @param aN Evaluated geometry interpolation
         */
        void evaluate_geometry_interpolation( const Matrix< DDRMat >& aXi, Matrix< DDRMat>& aN )
        {
            MORIS_ERROR( false, "Don't know how to evaluate a geometry interpolation for a T-matrix of dimension %u", N );
        }

        /**
         * Evaluates a shape function at a given point
         *
         * @param aXi Local coordinates
         * @param aN Evaluated shape function
         */
        void evaluate_shape_function( const Matrix< DDRMat >& aXi, Matrix< DDRMat >& aN )
        {
            MORIS_ERROR( false, "Don't know how to evaluate a shape function for a T-matrix of dimension %u", N );
        }

        /**
         * Gets the local coordinates of the corner nodes for a given child
         *
         * @param aChildIndex Child index
         * @param aXi Local coordinates
         */
        void get_child_corner_nodes( uint aChildIndex, Matrix< DDRMat >& aXi )
        {
            MORIS_ERROR( false, "Don't know how to get the corner nodes for a T-matrix of dimension %u", N );
        }

        /**
         * Gets the supporting points on another mesh of a given order
         *
         * @param aOrder Polynomial order
         * @return Supporting points
         */
        Matrix< DDRMat > get_supporting_points( uint aOrder )
        {
            MORIS_ERROR( false, "Don't know how to get supporting points for a T-matrix of dimension %u", N );
            return {0, 0};
        }

        /**
         * Evaluates the 2D/3D truncation weights based on the 1D weights and stores them internally
         *
         * @param aWeights 1D truncation weights
         */
        void evaluate_truncation_weights( const Matrix< DDRMat >& aWeights )
        {
            MORIS_ERROR( false, "Don't know how to evaluate truncation weights for a T-matrix of dimension %u", N );
        }

        /**
         * Populates the child matrices based on the given left and right matrix
         *
         * @param aTL Left matrix
         * @param aTR Right matrix
         */
        void populate_child_matrices( const Matrix< DDRMat >& aTL, const Matrix< DDRMat >& aTR )
        {
            MORIS_ERROR( false, "Don't know how to populate child matrices for a T-matrix of dimension %u", N );
        }

    };

    /**
     * Declare template specializations for 2D and 3D
     */
    template<> void T_Matrix< 2 >::evaluate_geometry_interpolation( const Matrix< DDRMat >& aXi, Matrix< DDRMat>& aN );
    template<> void T_Matrix< 3 >::evaluate_geometry_interpolation( const Matrix< DDRMat >& aXi, Matrix< DDRMat>& aN );
    template<> void T_Matrix< 2 >::evaluate_shape_function( const Matrix< DDRMat >& aXi, Matrix< DDRMat >& aN );
    template<> void T_Matrix< 3 >::evaluate_shape_function( const Matrix< DDRMat >& aXi, Matrix< DDRMat >& aN );
    template<> void T_Matrix< 2 >::get_child_corner_nodes( uint aChildIndex, Matrix< DDRMat >& aXi );
    template<> void T_Matrix< 3 >::get_child_corner_nodes( uint aChildIndex, Matrix< DDRMat >& aXi );
    template<> Matrix< DDRMat > T_Matrix< 2 >::get_supporting_points( uint aOrder );
    template<> Matrix< DDRMat > T_Matrix< 3 >::get_supporting_points( uint aOrder );
    template<> void T_Matrix< 2 >::evaluate_truncation_weights( const Matrix< DDRMat>& aWeights );
    template<> void T_Matrix< 3 >::evaluate_truncation_weights( const Matrix< DDRMat>& aWeights );
    template<> void T_Matrix< 2 >::populate_child_matrices( const Matrix< DDRMat >& aTL, const Matrix< DDRMat>& aTR );
    template<> void T_Matrix< 3 >::populate_child_matrices( const Matrix< DDRMat >& aTL, const Matrix< DDRMat>& aTR );

}

#endif /* SRC_HMR_CL_HMR_T_MATRIX_HPP_ */

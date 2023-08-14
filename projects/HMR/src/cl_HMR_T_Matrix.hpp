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
#include "fn_HMR_calculate_basis_identifier.hpp"
#include "cl_Matrix.hpp" //LINALG/src
#include "fn_eye.hpp"
#include "fn_HMR_bspline_shape.hpp"

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
                this->init_unity_matrix();
                this->init_basis_index();
                this->init_lagrange_matrix();
                if ( aTruncate )
                {
                    aBSplineMesh->evaluate_truncation_weights( mTruncationWeights );
                }
                else
                {
                    mBSplineMesh->evaluate_child_matrices( mBasisIndices, mChildMatrices );
                }
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
                const Element*                                    aRootBsplineElement,
                const Element*                                    aExtendedBsplineElement,
                moris::Cell< moris::Cell< const mtk::Vertex* > >& aRootBsplineBasis,
                moris::Cell< const mtk::Vertex* >&                aExtendedBsplineBasis,
                moris::Cell< Matrix< DDRMat > >&            aWeights )
        {
            // Background_Element_Base* aRootBSpBackgroundElement = aRootBsplineElement->get_background_element();
            // Background_Element_Base* aExtendedBSpBackgroundElement = aExtendedBsplineElement->get_background_element();

            const luint* tRootIJK     = aRootBsplineElement->get_ijk();
            const luint* tExtendedIJK = aExtendedBsplineElement->get_ijk();

            // number of bspline coefficients per direction
            uint tNodesPerDirection = mBSplineMesh->Mesh_Base::get_order() + 1;

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
            uint tNumberOfBases = mBSplineMesh->get_number_of_bases_per_element();

            // initialize the projection matrix with the correct size
            Matrix< DDRMat > tL2ProjectionMatrix( tNumberOfBases, tNumberOfBases );

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
                                tL2ProjectionMatrix( mBasisIndices( a ), mBasisIndices( b ) ) = tMatrices1D( 0 )( i, k ) * tMatrices1D( 1 )( j, l );
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
                                        tL2ProjectionMatrix( mBasisIndices( a ), mBasisIndices( b ) ) = tMatrices1D( 0 )( i, l ) * tMatrices1D( 1 )( j, q ) * tMatrices1D( 2 )( k, p );
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
            moris::Cell< const Basis* > tRootBasis;
            moris::Cell< const Basis* > tExtendedBasis;

            //reserve enough memory for each of them
            tRootBasis.reserve(tNumberOfBases);
            tExtendedBasis.reserve(tNumberOfBases);

            // fill out the cell data
            for ( uint i = 0; i < tNumberOfBases; i++ )
            {
                // get the basis
                const Basis* tBasis = aRootBsplineElement->get_basis( i );

                //if it is active add it to the cell
                if ( tBasis->is_active() )
                {
                    tRootBasis.push_back( tBasis );
                }

                // get the basis
                const Basis* tBasisExtended = aExtendedBsplineElement->get_basis( i );

                //if it is active add it to the cell
                if ( tBasisExtended->is_active() )
                {
                    tExtendedBasis.push_back( tBasisExtended);
                }
            }

            // resize the output cells
            aWeights.resize(tNumberOfBases );
            aRootBsplineBasis.resize(tNumberOfBases );
            aExtendedBsplineBasis.resize(tNumberOfBases );

            // loop over the basis and eliminate the zero values
            for ( uint iExtendedBasisIndex = 0; iExtendedBasisIndex < tNumberOfBases; iExtendedBasisIndex++ )
            {
                // set size for each of the extended basis
                aExtendedBsplineBasis(iExtendedBasisIndex ) = tExtendedBasis(iExtendedBasisIndex );
                aWeights(iExtendedBasisIndex ).set_size(tNumberOfBases, 1 );
                aRootBsplineBasis(iExtendedBasisIndex ).resize(tNumberOfBases );

                // find all the basis and weights of the root that hve non-zero values
                uint tNonzeroCount = 0;
                for ( uint iRootBasisIndex = 0; iRootBasisIndex < tNumberOfBases; iRootBasisIndex++ )
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

    protected:
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

            // calculate number of basis per element
            uint tNumberOfBases = mBSplineMesh->get_number_of_bases_per_element();

            // initialize index matrix
            mBasisIndices.set_size( tNumberOfBases, 1 );
            mBSplineIJK.set_size( N, tNumberOfBases );

            // Array of ijk position of basis
            luint tIJK[ N ];
            luint tOffset[ N ];
            for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < tNumberOfBases; iBasisIndexInElement++ )
            {
                // Get IJK position from element
                tElement->get_ijk_of_basis( iBasisIndexInElement, tIJK );

                // Loop over dimensions
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    // Store IJK
                    mBSplineIJK( iDimension, iBasisIndexInElement ) = tIJK[ iDimension ];

                    // Get offset
                    tOffset[ iDimension ] = mBSplineMesh->get_order( iDimension ) + 1;
                }

                // Get basis index for storage in T-matrix
                luint tBasisIndex = calculate_basis_identifier< N >( tIJK, tOffset );

                // Assign basis index
                mBasisIndices( tBasisIndex ) = iBasisIndexInElement;
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
            uint tNumberOfBases = mBSplineMesh->get_number_of_bases_per_element();

            // Set unity and zero matrices
            eye( tNumberOfBases, tNumberOfBases, mEye );
            mZero.set_size( tNumberOfBases, 1, 0.0 );
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
                    mLagrangeParamModified( iDimension, iNodeIndex ) = 2 * tScale * tIJKValue - 1.0;
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
            uint tNumberOfBases = mBSplineMesh->get_number_of_bases_per_element();

            // get number of Lagrange nodes
            uint tNumberOfNodes = mLagrangeParam.n_cols();

            // initialize T-Matrix for B-Spline to Lagrange conversion
            mTMatrixLagrange.set_size( tNumberOfNodes, tNumberOfBases, 1 );

            // Loop over all dimensions
            for ( uint iDimensionIndex = 0; iDimensionIndex < N; iDimensionIndex++ )
            {
                // get order of B-Spline mesh
                uint tBSplineOrder = mBSplineMesh->get_order( iDimensionIndex );

                // loop over all B-Spline Basis
                for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBases; iBasisIndex++ )
                {
                    // Loop over all Lagrange nodes
                    for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; iNodeIndex++ )
                    {
                        mTMatrixLagrange( iNodeIndex, iBasisIndex ) *= bspline_shape(
                                tBSplineOrder,
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
                this->get_child_corner_nodes( iChildIndex, tCorners );

                for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; iNodeIndex++ )
                {
                    // evaluate shape function for "geometry"
                    this->evaluate_geometry_interpolation( mLagrangeParam.get_column( iNodeIndex ), tNGeo );

                    // get parameter coordinates
                    Matrix< DDRMat > tXi = tNGeo * tCorners;

                    // evaluate shape function
                    this->evaluate_shape_function( tXi, tN );

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
                    // copy coordinates into point
                    tXi.get_column( iNodeOtherMesh, tPoint );

                    // evaluate shape function
                    this->evaluate_shape_function( tPoint,tN );

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
            uint tNumberOfBases = mBSplineMesh->get_number_of_bases_per_element();

            // get number of Lagrange nodes
            uint tNumberOfNodes = mLagrangeParam.n_cols();

            // initialize T-Matrix for B-Spline to Lagrange conversion
            mTMatrixLagrangeModified.set_size( tNumberOfNodes, tNumberOfBases, 1 );

            // Loop over all dimensions
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                // get order of B-Spline mesh
                uint tBSplineOrder = mBSplineMesh->get_order( iDimension );

                // loop over all B-Spline Basis
                for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBases; iBasisIndex++ )
                {
                    // Loop over all Lagrange nodes
                    for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; iNodeIndex++ )
                    {
                        mTMatrixLagrangeModified( iNodeIndex, iBasisIndex ) *= bspline_shape_extended(
                                tBSplineOrder,
                                mBSplineIJK( iDimension, iBasisIndex ),
                                mLagrangeParamModified( iDimension, iNodeIndex ) );
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------

        void get_extension_matrix_1d( real const & aShift, Matrix< DDRMat >& aExtensionMatrix )
        {
            switch ( mBSplineMesh->Mesh_Base::get_order() )
            {
                case 1:
                    aExtensionMatrix = { { 1.0 - aShift, aShift }, { -aShift, 1.0 + aShift } };
                    break;
                case 2:
                    aExtensionMatrix = { { 0.5 * ( aShift - 2.0 ) * ( aShift - 1.0 ), aShift * ( -aShift + 2.0 ), 0.5 * aShift * ( aShift - 1.0 ) },    //
                            { 0.5 * aShift * ( aShift - 1.0 ), -( aShift - 1.0 ) * ( aShift + 1.0 ), 0.5 * aShift * ( aShift + 1 ) },                       //
                            { 0.5 * aShift * ( aShift + 1.0 ), -aShift * ( aShift + 2.0 ), 0.5 * ( aShift + 1 ) * ( aShift + 2 ) } };
                    break;
                case 3:
                    aExtensionMatrix = {
                        {                                                                                                             //
                                -1.0 / 6.0 * std::pow( aShift, 3.0 ) + 1.0 * std::pow( aShift, 2.0 ) - 11.0 / 6.0 * aShift + 1.0,     //
                                +0.5 * std::pow( aShift, 3.0 ) - 2.5 * std::pow( aShift, 2.0 ) + 3.0 * aShift + 0.0,                  //
                                -0.5 * std::pow( aShift, 3.0 ) + 2.0 * std::pow( aShift, 2.0 ) - 1.5 * aShift + 0.0,                  //
                                +1.0 / 6.0 * std::pow( aShift, 3.0 ) - 0.5 * std::pow( aShift, 2.0 ) + 1.0 / 3.0 * aShift + 0.0 },    //
                        {                                                                                                             //
                                -1.0 / 6.0 * std::pow( aShift, 3.0 ) + 0.5 * std::pow( aShift, 2.0 ) - 1.0 / 3.0 * aShift + 0.0,      //
                                +0.5 * std::pow( aShift, 3.0 ) - 1.0 * std::pow( aShift, 2.0 ) - 0.5 * aShift + 1.0,                  //
                                -0.5 * std::pow( aShift, 3.0 ) + 0.5 * std::pow( aShift, 2.0 ) + 1.0 * aShift + 0.0,                  //
                                +1.0 / 6.0 * std::pow( aShift, 3.0 ) - 0.0 * std::pow( aShift, 2.0 ) - 1.0 / 6.0 * aShift + 0.0 },    //
                        {                                                                                                             //
                                -1.0 / 6.0 * std::pow( aShift, 3.0 ) - 0.0 * std::pow( aShift, 2.0 ) + 1.0 / 6.0 * aShift + 0.0,      //
                                +0.5 * std::pow( aShift, 3.0 ) + 0.5 * std::pow( aShift, 2.0 ) - 1.0 * aShift + 0.0,                  //
                                -0.5 * std::pow( aShift, 3.0 ) - 1.0 * std::pow( aShift, 2.0 ) + 0.5 * aShift + 1.0,                  //
                                +1.0 / 6.0 * std::pow( aShift, 3.0 ) + 0.5 * std::pow( aShift, 2.0 ) + 1.0 / 3.0 * aShift + 0.0 },
                        {                                                                                                             //
                                -1.0 / 6.0 * std::pow( aShift, 3.0 ) - 0.5 * std::pow( aShift, 2.0 ) - 1.0 / 3.0 * aShift + 0.0,      //
                                +0.5 * std::pow( aShift, 3.0 ) + 2.0 * std::pow( aShift, 2.0 ) + 1.5 * aShift + 0.0,                  //
                                -0.5 * std::pow( aShift, 3.0 ) - 2.5 * std::pow( aShift, 2.0 ) - 3.0 * aShift + 0.0,                  //
                                +1.0 / 6.0 * std::pow( aShift, 3.0 ) + 1.0 * std::pow( aShift, 2.0 ) + 11.0 / 6.0 * aShift + 1.0 }
                        //
                    };
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

}

#endif /* SRC_HMR_CL_HMR_T_MATRIX_HPP_ */

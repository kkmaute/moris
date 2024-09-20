/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Mesh.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_MESH_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_MESH_HPP_

#include "cl_HMR_Background_Element_Base.hpp"    //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp"       //HMR/src
#include "cl_HMR_BSpline_Element.hpp"            //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"          //HMR/src
#include "cl_HMR_Parameters.hpp"                 //HMR/src
#include "fn_HMR_calculate_basis_identifier.hpp"
#include "HMR_Globals.hpp"       //HMR/src
#include "moris_typedefs.hpp"    //COR/src
#include "cl_Stopwatch.hpp"      //CHR/src

namespace moris::hmr
{
    /**
     * B-spline element class
     *
     * @param P Polynomial degree in x-direction
     * @param Q Polynomial degree in y-direction
     * @param R Polynomial degree in z-direction
     */
    template< uint P, uint Q, uint R >
    class BSpline_Mesh : public BSpline_Mesh_Base
    {
        //! Number of dimensions
        static constexpr uint N = ( P > 0 ) + ( Q > 0 ) + ( R > 0 );

        //! Number of bases per element
        static constexpr uint B = ( P + 1 ) * ( Q + 1 ) * ( R + 1 );

        //! Number of children per basis
        static constexpr uint C = ( P + 2 - ( P == 0 ) ) * ( Q + 2 - ( Q == 0 ) ) * ( R + 2 - ( R == 0 ) );

        //! Container of degrees, used to avoid repeated conditionals
        static constexpr uint PQR[ 3 ] = { P, Q, R };

        //! Lookup table containing offset for node IDs
        luint mBasisLevelOffset[ gMaxNumberOfLevels ];

        //! Lookup table containing number of elements per dimension for each level
        luint mNumberOfElementsPerDimensionIncludingAura[ gMaxNumberOfLevels ][ N ];

        //! Lookup table for node IDs
        luint mMySubdomainOffset[ gMaxNumberOfLevels ][ N ];

        //! Lookup table containing number of basis per dimension for each level
        luint mNumberOfBasisPerDimensionIncludingPadding[ gMaxNumberOfLevels ][ N ];

        //! number of basis on coarsest level
        luint mNumberOfCoarsestBasisOnProc[ N ] = { 0 };

      public:
        // ----------------------------------------------------------------------------

        /**
         * Constructor for Lagrange Mesh
         *
         * @param[in] aParameters       ref to container of user defined settings
         * @param[in] aBackgroundMesh pointer to background mesh
         *
         */
        BSpline_Mesh(
                const Parameters*     aParameters,
                Background_Mesh_Base* aBackgroundMesh,
                uint                  aActivationPattern,
                uint                  aMeshIndex )
                : BSpline_Mesh_Base(    // <-- this initializer only copies what is passed here into member data and does not process anything
                          aParameters,
                          aBackgroundMesh,
                          this->get_min_order(),
                          aActivationPattern,
                          B )
        {
            // trace this operation
            Tracer tTracer( "HMR", "B-Spline Mesh #" + std::to_string( aMeshIndex ), "Create" );
            MORIS_LOG_INFO(
                    "Creating B-spline mesh index #%i with polynomial order p=%i on pattern #%i",
                    aMeshIndex,
                    mOrder,
                    aActivationPattern );

            // set and store the mesh index
            this->set_index( aMeshIndex );

            // Calculate child stencil (for multigrid)
            this->calculate_child_stencil();

            // ask background mesh for number of elements per ijk-direction for each refinement level
            this->get_number_of_elements_per_dimension();

            // calculate lookup table mBasisLevelOffset (which indicates which index ranges may be occupied for BFs on different refinement levels)
            this->calculate_lookup_tables();

            // calculate subdomain offset
            this->calculate_subdomain_offset();

            // calculate any value that can change after refinement
            this->update_mesh();
        }

        // ----------------------------------------------------------------------------

        /**
         * Default destructor.
         */
        ~BSpline_Mesh() override
        {
            mActiveBasisOnProc.clear();
            this->delete_pointers();
        }

        // ----------------------------------------------------------------------------

        uint
        get_order( uint aDimensionIndex ) override
        {
            return PQR[ aDimensionIndex ];
        }

        // ----------------------------------------------------------------------------

        uint
        get_min_order() override
        {
            // Ignores zero
            return std::min( { P - 1, Q - 1, R - 1 } ) + 1;
        }

        // ----------------------------------------------------------------------------

        uint
        get_max_order() override
        {
            return std::max( { P, Q, R } );
        }

        // ----------------------------------------------------------------------------

        uint
        get_number_of_bases_per_element() override
        {
            return B;
        }

        // ----------------------------------------------------------------------------

        luint
        calculate_basis_id(
                uint         aLevel,
                const luint* aIJK ) override
        {
            if ( aLevel < gMaxNumberOfLevels )
            {
                luint tIJK[ N ];
                luint tDimensionOffset[ N ];
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tIJK[ iDimension ]             = aIJK[ iDimension ] + mMySubdomainOffset[ aLevel ][ iDimension ];
                    tDimensionOffset[ iDimension ] = mNumberOfBasisPerDimensionIncludingPadding[ aLevel ][ iDimension ];
                }
                return calculate_basis_identifier< N >( tIJK, tDimensionOffset ) + mBasisLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

        // ----------------------------------------------------------------------------

        void
        evaluate_child_matrices(
                const Matrix< DDUMat >&     aBasisIndices,
                Vector< Matrix< DDRMat > >& aChildMatrices ) override
        {
            // Number of children per element
            uint tNumberOfChildren = std::pow( 2, N );

            // Calculate weight scaling factor
            real tScale = 1.0;
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                tScale /= std::pow( 2, PQR[ iDimension ] );
            }

            // Loop over dimensions
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                // Number of B-spline coefficients in this direction
                uint tNumCoefficients = PQR[ iDimension ] + 1;

                // Create temporary factors matrix to help with left/right
                Matrix< DDRMat > tFactors( tNumCoefficients, tNumCoefficients + 1, 0.0 );

                // Set factors matrix
                for ( uint iCoefficient = 0; iCoefficient <= tNumCoefficients; iCoefficient++ )
                {
                    for ( uint iOrder = 0; iOrder <= PQR[ iDimension ]; iOrder++ )
                    {
                        uint k = PQR[ iDimension ] - 2 * iOrder + iCoefficient;
                        if ( k <= tNumCoefficients )
                        {
                            tFactors( iOrder, iCoefficient ) = nchoosek( tNumCoefficients, k );
                        }
                    }
                }

                // left and right matrices
                Matrix< DDRMat > TL( tNumCoefficients, tNumCoefficients, 0.0 );
                Matrix< DDRMat > TR( tNumCoefficients, tNumCoefficients, 0.0 );

                // Fill matrices
                for ( uint iOrder = 0; iOrder <= PQR[ iDimension ]; iOrder++ )
                {
                    TL.set_column( iOrder, tFactors.get_column( iOrder ) );
                    TR.set_column( iOrder, tFactors.get_column( iOrder + 1 ) );
                }

                // determine number of basis per element
                uint tNumberOfBases = this->get_number_of_bases_per_element();

                // empty matrix
                Matrix< DDRMat > tScaleMat( tNumberOfBases, tNumberOfBases, tScale );

                // container for child relation matrices ( transposed! )
                aChildMatrices.resize( tNumberOfChildren, tScaleMat );

                // Left and right matrices
                Vector< Matrix< DDRMat > > tT( 2, Matrix< DDRMat >( tNumCoefficients, tNumCoefficients ) );

                // Fill matrices
                for ( uint iCoefficient = 0; iCoefficient < tNumCoefficients; iCoefficient++ )
                {
                    tT( 0 ).set_column( iCoefficient, tFactors.get_column( iCoefficient ) );
                    tT( 1 ).set_column( iCoefficient, tFactors.get_column( iCoefficient + 1 ) );
                }

                // Tensor product
                for ( uint tChildCol = 0; tChildCol < B; tChildCol++ )
                {
                    for ( uint tChildRow = 0; tChildRow < B; tChildRow++ )
                    {
                        uint tTRow = ( tChildCol / (uint)std::pow( tNumCoefficients, iDimension ) ) % tNumCoefficients;
                        uint tTCol = ( tChildRow / (uint)std::pow( tNumCoefficients, iDimension ) ) % tNumCoefficients;
                        for ( uint iChildIndex = 0; iChildIndex < tNumberOfChildren; iChildIndex++ )
                        {
                            aChildMatrices( iChildIndex )( aBasisIndices( tChildRow ), aBasisIndices( tChildCol ) ) *= tT( ( iChildIndex / (uint)std::pow( 2, iDimension ) ) % 2 )( tTRow, tTCol );
                        }
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------

        void
        evaluate_truncation_weights( Matrix< DDRMat >& aTruncationWeights ) override
        {
            // Calculate scale factor
            real tScale = 1.0;
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                tScale /= std::pow( 2, PQR[ iDimension ] );
            }

            // Allocate weights with scale factor
            aTruncationWeights.set_size( C, 1, tScale );

            // Perform multiplication
            uint tTruncationWeightIndex = 0;
            for ( uint iChildI = 0; iChildI < P + 2 - ( P == 0 ); iChildI++ )
            {
                for ( uint iChildJ = 0; iChildJ < Q + 2 - ( Q == 0 ); iChildJ++ )
                {
                    for ( uint iChildK = 0; iChildK < R + 2 - ( R == 0 ); iChildK++ )
                    {
                        aTruncationWeights( tTruncationWeightIndex++ ) *=
                                nchoosek( P + 1, iChildI )
                                * nchoosek( Q + 1, iChildJ )
                                * nchoosek( R + 1, iChildK );
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------

      private:
        // ----------------------------------------------------------------------------

        /**
         *  Private function, creates the mNodeLevelOffset lookup table.
         *
         *  @return void
         */
        void
        calculate_lookup_tables()
        {
            // Number of basis on level 0
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                mNumberOfBasisPerDimensionIncludingPadding[ 0 ][ iDimension ] = mNumberOfElementsPerDimensionIncludingAura[ 0 ][ iDimension ] + PQR[ iDimension ];
            }

            // Basis level offset on level 0
            mBasisLevelOffset[ 0 ] = 0;

            for ( uint iLevel = 1; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                // Start counting number of bases
                luint tNumberOfBasis = 1;

                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    // Calculate number of basis on higher levels
                    mNumberOfBasisPerDimensionIncludingPadding[ iLevel ][ iDimension ] =
                            2 * mNumberOfBasisPerDimensionIncludingPadding[ iLevel - 1 ][ iDimension ] + PQR[ iDimension ];

                    // Calculate number of nodes on this level
                    tNumberOfBasis *= mNumberOfBasisPerDimensionIncludingPadding[ iLevel - 1 ][ iDimension ];
                }

                // Add number of nodes to offset table
                mBasisLevelOffset[ iLevel ] = mBasisLevelOffset[ iLevel - 1 ] + tNumberOfBasis;
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Private function calculates the mMySubdomainOffset lookup table
         *
         * @return void
         */
        void
        calculate_subdomain_offset()
        {
            Matrix< DDLUMat > tIJK = mBackgroundMesh->get_subdomain_offset_of_proc();

            for ( uint iLevel = 0; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    mMySubdomainOffset[ iLevel ][ iDimension ] = tIJK( iDimension, iLevel );
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Internal function. Asks the background mesh for number of elements
         * per direction, stores result in  mAuraNumberOfElementsPerDimension
         *
         * @return void
         */
        void
        get_number_of_elements_per_dimension()
        {
            // get elements per level from background mesh
            Matrix< DDLUMat > tMat = mBackgroundMesh->get_number_of_elements_per_direction();

            // convert matrix to fixed size array
            for ( uint iLevel = 0; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    mNumberOfElementsPerDimensionIncludingAura[ iLevel ][ iDimension ] = tMat( iDimension, iLevel );
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * calculates XZY coordinates for each basis
         *
         * @return void
         */
        void
        calculate_basis_coordinates() override
        {
            // report on this operation
            MORIS_LOG_INFO( "B-Spline Mesh #%i: Computing basis function coordinates", this->get_index() );

            // get domain dimensions from settings
            const Vector< real >& tDomainDimensions = mParameters->get_domain_dimensions();

            // get number of elements on coarsest level from settings
            const Vector< luint >& tNumberOfElements = mParameters->get_number_of_elements_per_dimension();

            // calculate step width
            real tDeltaX[ gMaxNumberOfLevels ][ N ];

            // calculate width for first level
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                tDeltaX[ 0 ][ iDimension ] = tDomainDimensions( iDimension ) / ( (real)( tNumberOfElements( iDimension ) ) );
            }

            // loop over all higher levels
            for ( uint iLevel = 1; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tDeltaX[ iLevel ][ iDimension ] = 0.5 * tDeltaX[ iLevel - 1 ][ iDimension ];
                }
            }

            // get domain offset
            Matrix< DDRMat > tParametersOffset = mParameters->get_domain_offset();

            // domain offset
            real tOffset[ N ];

            // get coordinates from background mesh
            Matrix< DDRMat > tOffsetCoords = mBackgroundMesh->get_domain_offset();

            // unflatten coordinates to a normal array
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                tOffset[ iDimension ] = tOffsetCoords( iDimension );
            }

            // coordinate shift for B-Spline
            // 0.5 - 0.5*P
            real tShift[ gMaxNumberOfLevels ][ N ];
            for ( uint iLevel = 0; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tShift[ iLevel ][ iDimension ] = 0.5 * ( PQR[ iDimension ] + 1 ) - std::pow( 2, iLevel ) * PQR[ iDimension ];
                }
            }

            // loop over all nodes
            for ( auto tBasis : mAllBasisOnProc )
            {
                // get ijk position of node
                const luint* tIJK = tBasis->get_ijk();

                // get level of node
                luint tLevel = tBasis->get_level();

                // array containing coordinate
                real tXYZ[ N ];

                // loop over all dimensions
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tXYZ[ iDimension ] = ( tShift[ tLevel ][ iDimension ] + (real)( tIJK[ iDimension ] + mMySubdomainOffset[ tLevel ][ iDimension ] ) )
                                               * tDeltaX[ tLevel ][ iDimension ]
                                       + tOffset[ iDimension ];
                }

                // write XYZ coordinate into node
                tBasis->set_xyz( tXYZ );
            }
        }

        // ----------------------------------------------------------------------------

        void
        create_basis_on_level_zero() override
        {
            // Ask mesh for relevant ijk positions
            Matrix< DDLUMat > tNumElementsPerDirection = mBackgroundMesh->get_number_of_elements_per_direction_on_proc();

            // Unroll min and max i and j
            luint tTotalNumberOfCoarsestBases = 1;
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                mNumberOfCoarsestBasisOnProc[ iDimension ] = tNumElementsPerDirection( iDimension, 0 ) + PQR[ iDimension ];
                tTotalNumberOfCoarsestBases *= mNumberOfCoarsestBasisOnProc[ iDimension ];
            }

            // Size array
            mAllCoarsestBasisOnProc.resize( tTotalNumberOfCoarsestBases, nullptr );

            // Populate array
            luint tIJK[ N ];
            luint tBasisIndex = 0;
            populate_bases< N >( tIJK, tBasisIndex );
        }

        // ----------------------------------------------------------------------------

        /**
         * Populates the container of coarse basis pointers based on IJK positions
         *
         * @param D Number of dimensions left to process
         * @param aIJK IJK position to fill and use
         * @param aBasisIndex Index in the basis container to fill next
         */
        template< uint D, std::enable_if_t< ( D > 0 ) >* = nullptr >
        void
        populate_bases( luint* aIJK, luint& aBasisIndex )
        {
            // Loop over IJK
            for ( uint iBF = 0; iBF < mNumberOfCoarsestBasisOnProc[ D - 1 ]; iBF++ )
            {
                // Assign this IJK value
                aIJK[ D - 1 ] = iBF;

                // Go to next dimension
                populate_bases< D - 1 >( aIJK, aBasisIndex );
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * 0 specialization for populating coarse basis container, creates a new B-spline basis
         *
         * @param D Number of dimensions left to process (0)
         * @param aIJK IJK position
         * @param aBasisIndex Index in the basis container to fill
         */
        template< uint D, std::enable_if_t< ( D == 0 ) >* = nullptr >
        void
        populate_bases( luint* aIJK, luint& aBasisIndex )
        {
            // Create new basis
            mAllCoarsestBasisOnProc( aBasisIndex++ ) = new BSpline< P, Q, R >( aIJK, 0, gNoProcOwner );
        }

        //------------------------------------------------------------------------------

        void
        link_basis_to_elements_on_level_zero() override
        {
            // loop over all elements
            for ( auto tElement : mAllCoarsestElementsOnProc )
            {
                // init basis container
                tElement->init_basis_container();

                // loop over all basis of this element
                for ( uint iBasisIndex = 0; iBasisIndex < B; iBasisIndex++ )
                {
                    // Get IJK position of this basis
                    luint tIJK[ N ];
                    tElement->get_ijk_of_basis( iBasisIndex, tIJK );

                    // Get basis index
                    luint tCoarseBasisIndex = calculate_basis_identifier< N >( tIJK, mNumberOfCoarsestBasisOnProc );

                    // Insert point to basis into element
                    tElement->insert_basis( iBasisIndex, mAllCoarsestBasisOnProc( tCoarseBasisIndex ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        preprocess_bases_from_level(
                Vector< Element* >& aElements,
                Vector< Basis* >&   aBasisFunctions ) override
        {
            // reset flags for basis
            for ( Element* tElement : aElements )
            {
                // loop over all basis from this element
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                    if ( tBasis != nullptr )
                    {
                        // un-flag this basis
                        tBasis->unflag();

                        // deactivate this basis
                        tBasis->unuse();

                        // counts this element
                        tBasis->increment_element_counter();
                    }
                }
            }

            // initialize basis counter
            luint tBasisCount = 0;

            // get my rank
            moris_id tMyRank = par_rank();

            // ======================
            // go over all elements on the current level and mark the basis functions interpolating into them for usage and count them
            for ( Element* tElement : aElements )
            {
                // loop over all basis from this element
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                    if ( tBasis != nullptr )
                    {
                        // test if basis has been counted - make sure to not count basis functions twice
                        if ( !tBasis->is_flagged() )
                        {
                            // count this basis
                            ++tBasisCount;

                            // flag this basis
                            tBasis->flag();
                        }
                    }
                }

                // use basis if element is owned
                if ( tElement->get_owner() == tMyRank )
                {
                    // loop over all basis from this element
                    for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( iBasisIndexInElement );
                        if ( tBasis != nullptr )
                        {
                            tBasis->use();
                        }
                    }
                }
            }    // end for: each element on current level

            // assign memory for basis container
            aBasisFunctions.resize( tBasisCount, nullptr );

            // reset counter
            tBasisCount = 0;

            // ======================
            // go over all elements on the current level and initialize the neighbor containers of the BFs interpolating into the refined elements
            for ( Element* iElement : aElements )
            {
                // loop over all basis from this element and un-flag them again
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis* tBasis = iElement->get_basis( iBasisIndexInElement );

                    if ( tBasis != nullptr )
                    {
                        // if this basis has been processed already, reset the flag
                        if ( tBasis->is_flagged() )
                        {
                            // copy pointer to basis
                            aBasisFunctions( tBasisCount++ ) = tBasis;

                            // unflag this basis
                            tBasis->unflag();
                        }
                    }
                }

                // initialize basis neighbor container on basis functions interpolating into refined elements
                if ( iElement->is_refined() )
                {
                    // loop over all basis
                    for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                    {
                        // get pointer to basis
                        Basis* tBasis = iElement->get_basis( iBasisIndexInElement );

                        if ( tBasis != nullptr )
                        {
                            // assign memory of neighbor container for this basis
                            tBasis->init_neighbor_container();
                        }
                    }
                }
            }    // end for: each element on current level

            // ======================
            // link elements back to basis functions
            for ( auto iBF : aBasisFunctions )
            {
                // initialize element container
                iBF->init_element_container();
            }

            // loop over all elements
            for ( Element* tElement : aElements )
            {
                // loop over all basis from this element
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                    if ( tBasis != nullptr )
                    {
                        // insert this element into basis
                        tBasis->insert_element( tElement );
                    }
                }
            }

            // delete_unused_bases (nice feature, not sure if worth the effort)
            // this->delete_unused_bases( aLevel, aBackgroundElements, aBasisFunctions );              // FIXME Saves Memory

            // ======================
            // link basis functions with neighboring BFs
            for ( Element* tElement : aElements )
            {
                // calculate basis neighbors
                if ( tElement->is_refined() )
                {
                    // determine the BF's neighboring BFs
                    tElement->link_basis_with_neighbors( mAllElementsOnProc );
                }
            }

            // reset flag of all basis functions
            for ( auto iBasisFunction : aBasisFunctions )
            {
                // initialize element container
                iBasisFunction->unflag();
            }
        }

        //------------------------------------------------------------------------------

        void
        determine_basis_state( Vector< Basis* >& aBases ) override
        {
            // loop over all basis functions parsed into this function and flag them as (de-)activated and as (non-)refined
            for ( Basis* iBasisFunction : aBases )
            {
                // only process basis that are used by this proc
                if ( iBasisFunction->is_used() )
                {
                    // test number of elements in basis function's support
                    uint tNumberOfElements = iBasisFunction->get_element_counter();

                    // if the basis function's support extends beyond the padding (i.e. it is not fully supported on any level) it is irrelevant
                    if ( tNumberOfElements < B )
                    {
                        // mark this basis as de-activated
                        iBasisFunction->unset_active_flag();
                    }
                    else    // Basis function is fully supported within domain + padding
                    {
                        // check if any element in the BF's support are neither active nor refined which also indicates irrelevance
                        // TODO: for the truncation refactor the condition is that any element within the BF's support is active to be considered
                        bool tHasDeactivatedElement = false;
                        for ( uint iElementIndex = 0; iElementIndex < B; iElementIndex++ )
                        {
                            Element* tElement = iBasisFunction->get_element( iElementIndex );
                            if ( tElement->is_neither_active_nor_refined() )
                            {
                                tHasDeactivatedElement = true;
                                break;
                            }
                        }

                        // if the basis function is not fully supported by active or refined elements, deactivate it
                        if ( tHasDeactivatedElement )
                        {
                            iBasisFunction->unset_active_flag();
                        }
                        else
                        {
                            bool tIsActive = false;

                            // consider BF active if any of the elements in the basis function's support on its level are active
                            for ( uint iElementIndex = 0; iElementIndex < B; iElementIndex++ )
                            {
                                Element* tElement = iBasisFunction->get_element( iElementIndex );
                                if ( tElement->is_active() )
                                {
                                    tIsActive = true;

                                    // break loop
                                    break;
                                }
                            }

                            // the BF is supported by some active background element(s), therefore it remains active
                            if ( tIsActive )
                            {
                                // flag this basis as active
                                iBasisFunction->set_active_flag();
                            }
                            else    // the BF interpolates only into de-activated elements, hence it must be refined and fully replaced by finer BFs
                            {
                                // flag this basis as refined
                                iBasisFunction->set_refined_flag();
                            }
                        }    // end if: BF interpolates into de-activated element
                    }    // end if: BF is relevant and not outside the domain
                }    // end if: BF is used on processor
            }    // end for: all basis functions parsed into function
        }    // end function: BSpline_Mesh::determine_basis_state()

        //------------------------------------------------------------------------------

        void
        link_bases_to_parents() override
        {
            // ask background mesh for max number of levels
            uint tMaxLevel = mBackgroundMesh->get_max_level();

            // Cell containing children
            Vector< Basis* > tChildren;

            // loop over all levels but the last
            for ( uint iLevel = 0; iLevel < tMaxLevel; iLevel++ )
            {
                // make children from last step to parents
                Vector< Basis* > tBasis;
                this->collect_bases_from_level( iLevel, tBasis );

                // loop over all basis on this level
                for ( auto tParent : tBasis )
                {
                    // test if parent has children
                    if ( tParent->has_children() )
                    {
                        // loop over all children
                        for ( uint iChildIndex = 0; iChildIndex < C; iChildIndex++ )
                        {
                            Basis* tChild = tParent->get_child( iChildIndex );

                            // pointer may be null because we deleted unused basis
                            if ( tChild != nullptr )
                            {
                                // increment parent counter for child
                                tChild->increment_parent_counter();
                            }
                        }
                    }
                }

                // loop over all basis on this level
                for ( auto tParent : tBasis )
                {
                    // test if parent has children
                    if ( tParent->has_children() )
                    {
                        // loop over all children
                        for ( uint iChildIndex = 0; iChildIndex < C; iChildIndex++ )
                        {
                            Basis* tChild = tParent->get_child( iChildIndex );

                            // pointer may be null because we deleted unused basis
                            if ( tChild != nullptr )
                            {
                                // copy pointer of parent to child
                                tParent->get_child( iChildIndex )->insert_parent( tParent );
                            }
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        delete_unused_bases(
                uint                                aLevel,
                Vector< Background_Element_Base* >& aBackgroundElements,
                Vector< Basis* >&                   aBasis ) override
        {
            // start timer
            tic tTimer;

            // reset counter for debugging output
            luint tDeleteCount = 0;

            if ( aLevel > 0 )
            {
                // step 1: remove basis from parents

                // collect basis from upper level
                Vector< Basis* > tParents;

                this->collect_bases_from_level( aLevel - 1, tParents );

                // loop over all parents
                for ( auto tParent : tParents )
                {
                    // loop over all children of this parent
                    for ( uint iChildIndex = 0; iChildIndex < C; iChildIndex++ )
                    {
                        // get pointer to child
                        Basis* tChild = tParent->get_child( iChildIndex );

                        // test if child exists
                        if ( tChild != nullptr )
                        {
                            // test if basis is not used
                            if ( !tChild->is_used() )
                            {
                                // write null into child
                                tParent->insert_child( iChildIndex, nullptr );
                            }
                        }
                    }
                }

                // tidy up memory
                tParents.clear();

                // step 2: remove basis from basis neighbors
                for ( auto tBasis : aBasis )
                {
                    // loop over all neighbors
                    // ( number of neighbors per basis = mNumberOfNeighborsPerElement)
                    for ( uint k = 0; k < mNumberOfNeighborsPerElement; ++k )
                    {
                        // get pointer to neighbor
                        Basis* tNeighbor = tBasis->get_neighbor( k );

                        // test if neighbor exists
                        if ( tNeighbor != nullptr )
                        {
                            // test if neighbor is not used
                            if ( !tNeighbor->is_used() )
                            {
                                // write nullptr into neighbor
                                tBasis->insert_neighbor( k, nullptr );
                            }
                        }
                    }
                }

                // step 3: remove basis from elements
                for ( auto tBackElement : aBackgroundElements )
                {
                    // get pointer to element
                    Element* tElement = mAllElementsOnProc( tBackElement->get_memory_index() );

                    // loop over all basis
                    for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                        // test if basis exists
                        if ( tBasis != nullptr )
                        {
                            // test if basis is not used
                            if ( !tBasis->is_used() )
                            {
                                // write null into element
                                tElement->insert_basis( iBasisIndexInElement, nullptr );
                            }
                        }
                    }
                }

                // step 4: count active basis

                // init counter
                luint tBasisCount = 0;
                for ( auto tBasis : aBasis )
                {
                    // test if basis is used
                    if ( tBasis->is_used() )
                    {
                        // increment counter
                        ++tBasisCount;
                    }
                }

                // initialize output array
                Vector< Basis* > tBasisOut( tBasisCount, nullptr );

                // reset counter
                tBasisCount = 0;

                for ( auto tBasis : aBasis )
                {
                    // test if basis is used
                    if ( tBasis->is_used() )
                    {
                        // increment counter
                        tBasisOut( tBasisCount++ ) = tBasis;
                    }
                    else
                    {
                        delete tBasis;
                        ++tDeleteCount;
                    }
                }

                // remember number of basis for verbose output
                luint tNumberOfAllBasis = aBasis.size();

                // move output
                aBasis.clear();
                aBasis = tBasisOut;

                // stop timer
                real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

                // print output
                MORIS_LOG_INFO( "%s Deleted %lu unused basis of %lu total on level %u.",
                        proc_string().c_str(),
                        (long unsigned int)tDeleteCount,
                        (long unsigned int)tNumberOfAllBasis,
                        (unsigned int)aLevel );

                MORIS_LOG_INFO( "Deleting took %5.3f seconds.",
                        (double)tElapsedTime / 1000 );
            }
        }

        //------------------------------------------------------------------------------

        void
        collect_bases_from_level(
                uint              aLevel,
                Vector< Basis* >& aBasis ) override
        {
            Vector< Element* > tElements;

            this->collect_active_and_refined_elements_from_level( aLevel, tElements );

            // reset flags for basis
            for ( Element* tElement : tElements )
            {
                // loop over all basis from this element
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                    // test if basis exists
                    if ( tBasis != nullptr )
                    {
                        // unflag this basis
                        tBasis->unflag();
                    }
                }
            }

            // initialize basis counter
            luint tBasisCount = 0;

            // count basis on this level
            for ( Element* tElement : tElements )
            {
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                    // test if basis exists
                    if ( tBasis != nullptr )
                    {
                        // test if basis has been counted
                        if ( !tBasis->is_flagged() )
                        {
                            // count this basis
                            ++tBasisCount;

                            // flag this basis
                            tBasis->flag();
                        }
                    }
                }
            }

            // assign memory for basis container
            aBasis.resize( tBasisCount, nullptr );

            // reset basis counter
            tBasisCount = 0;

            // add basis to container
            for ( Element* tElement : tElements )
            {
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                    // test if basis exists
                    if ( tBasis != nullptr )
                    {
                        // test if basis has been added
                        if ( tBasis->is_flagged() )
                        {
                            // count this basis
                            aBasis( tBasisCount++ ) = tBasis;

                            // flag this basis
                            tBasis->unflag();
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        flag_refined_basis_of_owned_elements() override
        {
            if ( mParameters->use_multigrid() )
            {
                auto tMyRank = par_rank();

                // loop over all elements
                for ( Element* tElement : mAllElementsOnProc )
                {
                    // element must be neither padding or deactivated
                    if ( tElement->get_owner() == tMyRank )
                    {
                        if ( tElement->is_refined() and !tElement->is_padding() )
                        {
                            // loop over all basis of this element
                            for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                            {
                                // get pointer to basis
                                Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                                // flag basis if it is refined
                                if ( tBasis->is_refined() )
                                {
                                    tBasis->flag();
                                }
                            }
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        /**
         * Calculates child stencil for multigrid
         */
        void
        calculate_child_stencil()
        {
            // allocate matrix
            mChildStencil.set_size( C, 1 );

            // Start child counter
            uint tChildCounter = 0;

            // Switch over dimensions
            switch ( N )
            {
                case ( 2 ):
                {
                    for ( uint iChildI = 0; iChildI < P + 2; iChildI++ )
                    {
                        for ( uint iChildJ = 0; iChildJ < Q + 2; iChildJ++ )
                        {
                            mChildStencil( tChildCounter++ ) = nchoosek( P + 1, iChildI )
                                                             * nchoosek( Q + 1, iChildJ )
                                                             / std::pow( 2, P )
                                                             / std::pow( 2, Q );
                        }
                    }
                    break;
                }
                case ( 3 ):
                {
                    for ( uint iChildI = 0; iChildI < P + 2; iChildI++ )
                    {
                        for ( uint iChildJ = 0; iChildJ < Q + 2; iChildJ++ )
                        {
                            for ( uint iChildK = 0; iChildK < R + 2; iChildK++ )
                            {
                                mChildStencil( tChildCounter++ ) = nchoosek( P + 1, iChildI )
                                                                 * nchoosek( Q + 1, iChildJ )
                                                                 * nchoosek( R + 1, iChildK )
                                                                 / std::pow( 2, P )
                                                                 / std::pow( 2, Q )
                                                                 / std::pow( 2, R );
                            }
                        }
                    }
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Invalid dimension." );
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Creates a new B-spline element based on the provided background element, N, and P
         *
         * @param aBackgroundElement Background element
         * @return Created B-spline element
         */
        Element*
        create_element( Background_Element_Base* aBackgroundElement ) override
        {
            // Check for dimension less than or equal to 3 and degree less than or equal to 3
            MORIS_ERROR( P <= 3 and Q <= 3 and R <= 3, "hmr::BSpline_Mesh::create_element() - Don't know how to create B-Spline element." );

            // Create element
            return new BSpline_Element< P, Q, R >( aBackgroundElement, mActivationPattern );
        }

        // ----------------------------------------------------------------------------
    };
}    // namespace moris::hmr

#endif /* SRC_HMR_CL_HMR_BSPLINE_MESH_HPP_ */

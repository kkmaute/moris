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

#include "cl_HMR_Background_Element_Base.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Element.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "HMR_Globals.hpp" //HMR/src
#include "typedefs.hpp" //COR/src
#include "cl_Stopwatch.hpp" //CHR/src

namespace moris::hmr
{
    /**
     * B-spline element class
     *
     * @tparam P Polynomial degree in x-direction
     * @tparam Q Polynomial degree in y-direction
     * @tparam R Polynomial degree in z-direction
     */
    template< uint P, uint Q, uint R >
    class BSpline_Mesh : public BSpline_Mesh_Base
    {
        //! Number of dimensions
        static constexpr uint N = ( P > 0 ) + ( Q > 0 ) + ( R > 0 );
        
        //! Number of bases per element
        static constexpr uint B = ( P + 1 ) * ( Q + 1 ) * ( R + 1 );

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
                const Parameters      * aParameters,
                Background_Mesh_Base  * aBackgroundMesh,
                uint aActivationPattern )
        : BSpline_Mesh_Base(
                aParameters,
                aBackgroundMesh,
                this->get_min_order(),
                aActivationPattern,
                B )
        {
            // ask background mesh for number of elements per ijk-direction
            this->get_number_of_elements_per_dimension();

            // calculate lookup table mBasisLevelOffset
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

        uint get_order( uint aDimensionIndex ) override
        {
            return PQR[ aDimensionIndex ];
        }

        // ----------------------------------------------------------------------------

        uint get_min_order() override
        {
            // Ignores zero
            return std::min( { P - 1, Q - 1, R - 1 } ) + 1;
        }

        // ----------------------------------------------------------------------------

        uint get_max_order() override
        {
            return std::max( { P, Q, R } );
        }

        // ----------------------------------------------------------------------------

        uint get_number_of_bases_per_element() override
        {
            return ( P + 1 ) * ( Q + 1 ) * ( R + 1 );
        }

    private:

        // ----------------------------------------------------------------------------

        luint calculate_basis_id(
                uint         aLevel,
                const luint* aIJK ) override
        {
            if ( aLevel < gMaxNumberOfLevels )
            {
                luint tIJK[ N ];
                luint tDimensionOffset[ N ];
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tIJK[ iDimension ] = aIJK[ iDimension ] + mMySubdomainOffset[ aLevel ][ iDimension ];
                    tDimensionOffset[ iDimension ] = mNumberOfBasisPerDimensionIncludingPadding[ aLevel ][ iDimension ];
                }
                return this->calculate_basis_identifier( tIJK, tDimensionOffset ) + mBasisLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

        /**
         *  Private function, creates the mNodeLevelOffset lookup table.
         *
         *  @return void
         */
        void calculate_lookup_tables()
        {
            // Number of basis on level 0
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                mNumberOfBasisPerDimensionIncludingPadding[ 0 ][ iDimension ] = mNumberOfElementsPerDimensionIncludingAura[ 0 ][ iDimension ] + PQR[ iDimension ];
            }

            // Basis level offset on level 0
            mBasisLevelOffset[ 0 ] = 0;

            for( uint iLevel = 1; iLevel < gMaxNumberOfLevels; iLevel++ )
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
        void calculate_subdomain_offset()
        {
            Matrix< DDLUMat > tIJK = mBackgroundMesh->get_subdomain_offset_of_proc();

            for( uint iLevel = 0; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    mMySubdomainOffset[ iLevel ][ iDimension ] = tIJK(iDimension, iLevel );
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
        void get_number_of_elements_per_dimension()
        {
            // get elements per level from background mesh
            Matrix< DDLUMat > tMat = mBackgroundMesh->get_number_of_elements_per_direction();

            // convert matrix to fixed size array
            for( uint iLevel = 0; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    mNumberOfElementsPerDimensionIncludingAura[ iLevel ][ iDimension ] = tMat(iDimension, iLevel );
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * calculates XZY coordinates for each basis
         *
         * @return void
         */
        void calculate_basis_coordinates() override
        {
            // get domain dimensions from settings
            Matrix< DDRMat > tDomainDimensions = mParameters->get_domain_dimensions();

            // get number of elements on coarsest level from settings
            Matrix< DDLUMat > tNumberOfElements = mParameters->get_number_of_elements_per_dimension();

            // calculate step width
            real tDeltaX[ gMaxNumberOfLevels ][ N ];

            // calculate width for first level
            for( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                tDeltaX[ 0 ][ iDimension ] = tDomainDimensions(iDimension ) / ( ( real ) ( tNumberOfElements(iDimension ) ) );
            }

            // loop over all higher levels
            for( uint iLevel = 1; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tDeltaX[ iLevel ][ iDimension ] = 0.5 * tDeltaX[ iLevel-1 ][ iDimension ];
                }
            }

            // get domain offset
            Matrix< DDRMat > tParametersOffset = mParameters->get_domain_offset();

            // domain offset
            real tOffset[ N ];

            // get coordinates from background mesh
            Matrix< DDRMat > tOffsetCoords = mBackgroundMesh->get_domain_offset();

            // unflatten coordinates to a normal array
            for( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                tOffset[ iDimension ] = tOffsetCoords(iDimension );
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
            for( auto tBasis : mAllBasisOnProc )
            {
                // get ijk position of node
                const luint* tIJK = tBasis->get_ijk();

                // get level of node
                luint tLevel = tBasis->get_level();

                // array containing coordinate
                real tXYZ[ N ];

                // loop over all dimensions
                for( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tXYZ[ iDimension ] = ( tShift[ tLevel ][ iDimension ] + ( real ) ( tIJK[ iDimension ]
                                                                      + mMySubdomainOffset[ tLevel ][ iDimension ] ) )
                                                                      * tDeltaX[ tLevel ][ iDimension ] + tOffset[ iDimension ];
                }

                // write XYZ coordinate into node
                tBasis->set_xyz( tXYZ );
            }
        }

        // ----------------------------------------------------------------------------

        void create_basis_on_level_zero() override
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

        /**
         * Populates the container of coarse basis pointers based on IJK positions
         *
         * @tparam D Number of dimensions left to process
         * @param aIJK IJK position to fill and use
         * @param aBasisIndex Index in the basis container to fill next
         */
        template< uint D, std::enable_if_t< ( D > 0 ) >* = nullptr >
        void populate_bases( luint* aIJK, luint& aBasisIndex )
        {
            // Loop over IJK
            for ( uint i = 0; i < mNumberOfCoarsestBasisOnProc[ D - 1 ]; i++ )
            {
                // Assign this IJK value
                aIJK[ D - 1 ] = i;

                // Go to next dimension
                populate_bases< D - 1 >( aIJK, aBasisIndex );
            }
        }

        /**
         * 0 specialization for populating coarse basis container, creates a new B-spline basis
         *
         * @tparam D Number of dimensions left to process (0)
         * @param aIJK IJK position
         * @param aBasisIndex Index in the basis container to fill
         */
        template< uint D, std::enable_if_t< ( D == 0 ) >* = nullptr >
        void populate_bases( luint* aIJK, luint& aBasisIndex )
        {
            // Create new basis
            mAllCoarsestBasisOnProc( aBasisIndex++ ) = new BSpline< P, Q, R >( aIJK, 0, gNoProcOwner );
        }

        //------------------------------------------------------------------------------

        void link_basis_to_elements_on_level_zero() override
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
                    luint tCoarseBasisIndex = this->calculate_basis_identifier( tIJK, mNumberOfCoarsestBasisOnProc );

                    // Insert point to basis into element
                    tElement->insert_basis( iBasisIndex, mAllCoarsestBasisOnProc( tCoarseBasisIndex ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void preprocess_bases_from_level(
                Cell< Element * > & aElements,
                Cell< Basis * >   & aBasis ) override
        {
            // reset flags for basis
            for ( Element * tElement : aElements )
            {
                // loop over all basis from this element
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis * tBasis = tElement->get_basis( iBasisIndexInElement );

                    if ( tBasis != nullptr )
                    {
                        // unflag this basis
                        tBasis->unflag();

                        // unuse this basis
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

            for ( Element * tElement : aElements )
            {
                // loop over all basis from this element
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis * tBasis = tElement->get_basis( iBasisIndexInElement );

                    if ( tBasis != nullptr )
                    {
                        // test if basis has been counted
                        if ( ! tBasis->is_flagged() )
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
                        Basis * tBasis = tElement->get_basis( iBasisIndexInElement );
                        if ( tBasis != nullptr )
                        {
                            tBasis->use();
                        }
                    }
                }
            }

            // assign memory for basis container
            aBasis.resize( tBasisCount, nullptr );

            // reset counter
            tBasisCount = 0;

            // loop over all elements
            for ( Element * tElement : aElements )
            {
                // loop over all basis from this element
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis * tBasis = tElement->get_basis( iBasisIndexInElement );

                    if ( tBasis != nullptr )
                    {
                        // test if this basis has been processed already
                        if ( tBasis->is_flagged() )
                        {
                            // copy pointer to basis
                            aBasis( tBasisCount++ ) = tBasis;

                            // unflag this basis
                            tBasis->unflag();
                        }
                    }
                }

                // init neighbor container
                if ( tElement->is_refined() )
                {
                    // loop over all basis
                    for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                    {
                        // get pointer to basis
                        Basis * tBasis = tElement->get_basis( iBasisIndexInElement );

                        if ( tBasis != nullptr )
                        {
                            // assign memory of neighbor container for this basis
                            tBasis->init_neighbor_container();
                        }
                    }
                }
            } // end loop over elements

            // link elements to basis
            for ( auto tBasis : aBasis )
            {
                // initialize element container
                tBasis->init_element_container();
            }

            // loop over all elements
            for ( Element * tElement : aElements )
            {
                // loop over all basis from this element
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    Basis * tBasis = tElement->get_basis( iBasisIndexInElement );

                    if ( tBasis != nullptr )
                    {
                        // insert this element into basis
                        tBasis->insert_element( tElement );
                    }
                }
            }

            // delete_unused_bases (nice feature, not sure if worth the effort)
            //this->delete_unused_bases( aLevel, aBackgroundElements, aBasis );              // FIXME Saves Memory

            // link basis with neighbors
            for ( Element * tElement : aElements )
            {
                // calculate basis neighborship
                if ( tElement->is_refined() )
                {
                    // determine basis neighborship
                    tElement->link_basis_with_neighbors( mAllElementsOnProc );
                }
            }

            // reset flag of all basis
            for ( auto tBasis : aBasis )
            {
                // initialize element container
                tBasis->unflag();
            }
        }

        //------------------------------------------------------------------------------

        void delete_unused_bases(
                uint                               aLevel,
                Cell< Background_Element_Base* > & aBackgroundElements,
                Cell< Basis* >                   & aBasis ) override
        {
            // start timer
            tic tTimer;

            // reset counter for debugging output
            luint tDeleteCount = 0;

            if ( aLevel > 0 )
            {
                // step 1: remove basis from parents

                // collect basis from upper level
                Cell< Basis* > tParents;

                this->collect_bases_from_level( aLevel - 1, tParents );

                // loop over all parents
                for ( auto tParent : tParents )
                {
                    // loop over all children of this parent
                    for ( uint k=0; k<mNumberOfChildrenPerBasis; ++k )
                    {
                        // get pointer to child
                        Basis* tChild = tParent->get_child( k );

                        // test if child exists
                        if ( tChild != nullptr )
                        {
                            // test if basis is not used
                            if ( ! tChild->is_used() )
                            {
                                // write null into child
                                tParent->insert_child( k, nullptr );
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
                    for ( uint k=0; k<mNumberOfNeighborsPerElement; ++k )
                    {
                        // get pointer to neighbor
                        Basis* tNeighbor = tBasis->get_neighbor( k );

                        // test if neighbor exists
                        if ( tNeighbor != nullptr )
                        {
                            // test if neighbor is not used
                            if ( ! tNeighbor->is_used() )
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
                            if ( ! tBasis->is_used() )
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
                for ( auto tBasis: aBasis )
                {
                    // test if basis is used
                    if ( tBasis->is_used() )
                    {
                        // increment counter
                        ++tBasisCount;
                    }
                }

                // initialize output array
                Cell< Basis* > tBasisOut( tBasisCount, nullptr );

                // reset counter
                tBasisCount = 0;

                for ( auto tBasis: aBasis )
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
                aBasis = std::move( tBasisOut );

                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "%s Deleted %lu unused basis of %lu total on level %u.",
                                proc_string().c_str(),
                                ( long unsigned int ) tDeleteCount,
                                ( long unsigned int ) tNumberOfAllBasis,
                                ( unsigned int )      aLevel);

                MORIS_LOG_INFO( "Deleting took %5.3f seconds.",
                                ( double ) tElapsedTime / 1000 );
            }
        }

        //------------------------------------------------------------------------------

        void collect_bases_from_level(
                uint             aLevel,
                Cell< Basis* > & aBasis ) override
        {
            Cell< Element * > tElements;

            this->collect_active_and_refined_elements_from_level( aLevel, tElements );

            // reset flags for basis
            for ( Element * tElement : tElements )
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
            for ( Element * tElement : tElements )
            {
                for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( iBasisIndexInElement );

                    // test if basis exists
                    if ( tBasis != nullptr )
                    {
                        // test if basis has been counted
                        if ( ! tBasis->is_flagged() )
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
            for ( Element * tElement : tElements )
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

        void flag_refined_basis_of_owned_elements() override
        {
            if ( mParameters->use_multigrid() )
            {
                auto tMyRank = par_rank();

                // loop over all elements
                for ( Element * tElement : mAllElementsOnProc )
                {
                    // element must be neither padding or deactive
                    if ( tElement->get_owner() == tMyRank )
                    {
                        if ( tElement->is_refined() and ! tElement->is_padding() )
                        {
                            // loop over all basis of this element
                            for ( uint iBasisIndexInElement = 0; iBasisIndexInElement < B; iBasisIndexInElement++ )
                            {
                                // get pointer to basis
                                Basis * tBasis = tElement->get_basis( iBasisIndexInElement );

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

        // ----------------------------------------------------------------------------

        /**
         * Creates a new B-spline element based on the provided background element, N, and P
         *
         * @param aBackgroundElement Background element
         * @return Created B-spline element
         */
        Element * create_element( Background_Element_Base* aBackgroundElement ) override
        {
            // Check for dimension less than or equal to 3 and degree less than or equal to 3
            MORIS_ERROR( P <= 3 and Q <= 3 and R <= 3, "Don't know how to create B-Spline element." );

            // Create element
            return new BSpline_Element< P, Q, R >( aBackgroundElement, mActivationPattern );
        }

        // ----------------------------------------------------------------------------

        /**
         * Calculates a unique basis identifier for a given IJK position and offsets.
         * Can be an index or an ID depending on the offset.
         *
         * @tparam D Number of dimensions of the IJK and offset arrays
         * @param aIJK IJK position
         * @param aDimensionOffset Offset array for each dimension
         * @return Unique basis identifier
         */
        template< uint D = N >
        static luint calculate_basis_identifier(
                const luint* aIJK,
                const luint* aDimensionOffset )
        {
            luint tIdentifier = 0;
            for ( uint iDimension = 0; iDimension < D; iDimension++)
            {
                luint tOffsetTerm = aIJK[ iDimension ];
                for ( uint iPreviousDimension = 0; iPreviousDimension < iDimension; iPreviousDimension++ )
                {
                    tOffsetTerm *= aDimensionOffset[ iPreviousDimension ];
                }
                tIdentifier += tOffsetTerm;
            }
            return tIdentifier;
        }

        // ----------------------------------------------------------------------------

    };
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_MESH_HPP_ */


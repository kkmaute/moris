/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_HPP_

#include "cl_HMR_Basis.hpp"    //HMR/src

namespace moris::hmr
{

    /**
     * B-spline class
     *
     * @tparam P Polynomial degree in x-direction
     * @tparam Q Polynomial degree in y-direction
     * @tparam R Polynomial degree in z-direction
     */
    template< uint P, uint Q, uint R >
    class BSpline : public Basis
    {
        //! Number of dimensions
        static constexpr uint N = ( P > 0 ) + ( Q > 0 ) + ( R > 0 );

        //! Number of children
        static constexpr uint C = ( P + 2 - ( P == 0 ) ) * ( Q + 2 - ( Q == 0 ) ) * ( R + 2 - ( R == 0 ) );

        //! Number of neighbors
        static constexpr uint B = ( N == 1 ? 2 : ( N == 2 ? 8 : 26 ) );

        //! memory position in active cell
        luint mActiveIndex = gNoEntityID;

        //! local ijk position on proc
        luint mIJK[ N ];

#ifdef MORIS_HAVE_DEBUG
        //! global coordinates
        real mXYZ[ N ] = { 0 };
#endif

        //! flag telling if the basis is active
        bool mActiveFlag = false;

        //! flag telling if the basis is refined
        bool mRefinedFlag = false;

        //! flag that tells if the Neighbor array is allocated
        bool mNeighborsFlag = false;

        //! flag that tells of the children array is allocated
        bool mChildrenFlag = false;

        //! container with neighbors on the same level
        Basis** mNeighbors;
        // Basis* mNeighbors[ B ] = { nullptr };

        //! container with children
        Basis** mChildren;
        // Basis* mChildren [ C ] = { nullptr };

        //! Parent counter. Only needed for testing.
        uint mNumberOfParents = 0;
        uint mParentCounter   = 0;

        //! container for parents
        Basis** mParents;

        //! counter for connected Basis
        uint mNumberOfConnectedBasis = 0;

        //! flag telling if the basis container is allocated
        bool mConnectedFlag = false;

        //! container for connected basis
        Basis** mConnectedBasis;

        //------------------------------------------------------------------------------

      public:
        // -----------------------------------------------------------------------------

        /**
         * default constructor
         *
         * @param[in]   aIJK        ijk position of node
         * @param[in]   aLevel      level on which basis exists
         * @param[in]   aOwner      owner of basis
         */
        BSpline(
                const luint* aIJK,
                uint         aLevel,
                uint         aOwner )
                : Basis( aLevel, aOwner )
        {
            // save ijk position in memory.
            for ( uint k = 0; k < N; ++k )
            {
                mIJK[ k ] = aIJK[ k ];
            }

            // set used flag on by default
            mUsedFlag = true;
        }

        //------------------------------------------------------------------------------

        /**
         * Destructor
         */
        ~BSpline()
        {
            // test if basis has elements
            if ( mNumberOfConnectedElements != 0 )
            {
                mElements.clear();
            }

            // test if this basis has neighbors
            if ( mNeighborsFlag )
            {
                delete[] mNeighbors;
            }

            // test if this basis has children
            if ( mChildrenFlag )
            {
                delete[] mChildren;
            }

            // test if parent container is used
            if ( mNumberOfParents != 0 )
            {
                delete[] mParents;
            }

            if ( mConnectedFlag )
            {
                delete[] mConnectedBasis;
            }
        };

        //------------------------------------------------------------------------------

        /**
         * Sets the state of this basis to "active".
         *
         * @return void
         */
        void
        set_active_flag()
        {
            // set active flag on
            mActiveFlag = true;

            // an active element can not be refined at the same time
            mRefinedFlag = false;
        }

        //------------------------------------------------------------------------------

        /**
         *Sets the state of this basis to "refined".
         *
         * @return void
         */
        void
        set_refined_flag()
        {
            // a refined element is not active
            mActiveFlag = false;

            // set element as refined
            mRefinedFlag = true;
        }

        //------------------------------------------------------------------------------

        /**
         *Sets the state of this basis to "deactivated".
         *
         * @return void
         */
        void
        unset_active_flag()
        {
            // a refined element is not active
            mActiveFlag = false;

            // set element as refined
            mRefinedFlag = false;
        }

        //------------------------------------------------------------------------------

        /**
         * tells if a Basis is active
         *
         * @return bool   true if active
         */
        bool
        is_active() const
        {
            return mActiveFlag;
        }

        //------------------------------------------------------------------------------

        /**
         * tells if a Basis is refined
         *
         * @return bool   true if refined
         */
        bool
        is_refined()
        {
            return mRefinedFlag;
        }

        //------------------------------------------------------------------------------

        /**
         * Returns an array of size [N] telling the proc local ijk-position
         * of the B-Spline on the current level.
         *
         * @return luint pointer to array containing ijk-position
         *               careful: node must not go out of scope.
         */
        const luint*
        get_ijk() const
        {
            return mIJK;
        }

        // -----------------------------------------------------------------------------

        /**
         * set XYZ coordinates
         *
         * @param[in] aXYZ    array containing coordinates
         *
         * @return void
         */
        void
        set_xyz( const real* aXYZ )
        {
#ifdef MORIS_HAVE_DEBUG
            // save ijk position in memory.
            for ( uint k = 0; k < N; ++k )
            {
                mXYZ[ k ] = aXYZ[ k ];
            }
#endif
        }

        // -----------------------------------------------------------------------------

        /**
         * Returns an array of size [N] telling the xyz-position
         * of the B-Spline
         *
         * @return double pointer to array containing xyz-position
         *               careful: node must not go out of scope.
         */
        const real*
        get_xyz() const
        {
#ifdef MORIS_HAVE_DEBUG
            return mXYZ;
#else
            MORIS_ERROR( false, "get_xyz() If you wish this function to work for B-Splines and non DEBUG, delete the #ifdef MORIS_HAVE_DEBUG around mXYZ in the class BSplin" );
            return nullptr;
#endif
        }

        // -----------------------------------------------------------------------------

        /**
         * MTK Interface: return the coords of this node as Moris::Mat
         */
#ifdef MORIS_HAVE_DEBUG
        Matrix< DDRMat >
        get_coords() const
        {
            Matrix< DDRMat > aCoords( 1, N );
            for ( uint k = 0; k < N; ++k )
            {
                aCoords( k ) = mXYZ[ k ];
            }
            return aCoords;
        }
#endif

        // -----------------------------------------------------------------------------

        /**
         * reserves the memory for the  neighbor container
         *
         * @return void
         */
        void
        init_neighbor_container()
        {
            if ( !mNeighborsFlag )
            {
                // reserve array
                mNeighbors = new Basis*[ B ];

                // reset array
                for ( uint iBF = 0; iBF < B; ++iBF )
                {
                    mNeighbors[ iBF ] = nullptr;
                }

                // set flag
                mNeighborsFlag = true;
            }
        }

        // -----------------------------------------------------------------------------

        void
        delete_neighbor_container()
        {
            if ( mNeighborsFlag )
            {
                delete[] mNeighbors;
                mNeighborsFlag = false;
            }
        }

        // -----------------------------------------------------------------------------

        void
        init_connection_container()
        {
            if ( !mConnectedFlag )
            {
                mConnectedBasis = new Basis*[ mNumberOfConnectedBasis ];

                // reset array
                for ( uint k = 0; k < mNumberOfConnectedBasis; ++k )
                {
                    mConnectedBasis[ k ] = nullptr;
                }

                // reset the counter
                mNumberOfConnectedBasis = 0;

                // set the flag
                mConnectedFlag = true;
            }
        }

        // -----------------------------------------------------------------------------

        void
        delete_connection_container()
        {
            if ( mConnectedFlag )
            {
                delete[] mConnectedBasis;
                mConnectedFlag = false;
            }
        }

        // ----------------------------------------------------------------------------

        void
        increment_connection_counter()
        {
            ++mNumberOfConnectedBasis;
        }
        // -----------------------------------------------------------------------------

        void
        insert_connected_basis( Basis* aBasis )
        {
            mConnectedBasis[ mNumberOfConnectedBasis++ ] = aBasis;
        }

        // -----------------------------------------------------------------------------

        Basis*
        get_connected_basis( uint aBasisNumber )
        {
            return mConnectedBasis[ aBasisNumber ];
        }

        // -----------------------------------------------------------------------------

        const Basis*
        get_connected_basis( uint aBasisNumber ) const
        {
            return mConnectedBasis[ aBasisNumber ];
        }

        // -----------------------------------------------------------------------------

        uint
        get_number_of_connected_basis() const
        {
            return mNumberOfConnectedBasis;
        }

        // -----------------------------------------------------------------------------

        /**
         * Tells if children of this basis have been processed already
         */
        bool
        has_children()
        {
            return mChildrenFlag;
        }

        // -----------------------------------------------------------------------------

        /**
         * reserves the memory and initializes it with nullptr for the  neighbor container.
         * The children flag is set to true
         *
         * @return void
         */
        void
        init_children_container()
        {
            // reserve array
            mChildren = new Basis*[ C ];

            // reset array
            for ( uint k = 0; k < C; ++k )
            {
                mChildren[ k ] = nullptr;
            }

            // set flag
            mChildrenFlag = true;
        }

        uint
        get_number_of_children() override
        {
            return mChildrenFlag * C;
        }

        // -----------------------------------------------------------------------------

        /**
         * inserts a child to the specified position
         *
         * @param[ in ] aChildNumber    index of child
         * @param[ in ] aChild          pointer to child
         *
         * @return void
         */
        void
        insert_child( uint aChildNumber,
                Basis*     aChild )
        {
            mChildren[ aChildNumber ] = aChild;
        }

        // -----------------------------------------------------------------------------

        /**
         * returns the pointer to a child at specified position
         *
         * @param[ in ] aChildNumber    index of child
         *
         * @return Basis*               pointer to child
         */
        Basis*
        get_child( uint aChildNumber )
        {
            if ( mChildrenFlag )
            {
                return mChildren[ aChildNumber ];
            }
            else
            {
                return nullptr;
            }
        }

        // -----------------------------------------------------------------------------

        /**
         * Returns the specified neighbor of a basis.
         * Numbering scheme is the same as used for elements.
         * Returns nullptr if neighbor does not exist
         *
         * @param[ in ] aNeighborNumber    index of neighbor
         *
         * @return Basis*
         */
        Basis*
        get_neighbor( uint aNeighborNumber )
        {
            if ( mNeighborsFlag )
            {
                return mNeighbors[ aNeighborNumber ];
            }
            else
            {
                return nullptr;
            }
        }

        // -----------------------------------------------------------------------------

        void
        insert_neighbor( uint aNeighborNumber,
                Basis*        aNeighbor )
        {
            MORIS_ASSERT( mNeighborsFlag, "Can't insert neighbor if container is not set" );
            mNeighbors[ aNeighborNumber ] = aNeighbor;
        }

        // -----------------------------------------------------------------------------

        void
        flag_descendants()
        {
            if ( mFlag == false )
            {
                // flag myself
                mFlag = true;

                // flag children
                // test if children exist
                if ( mChildrenFlag )
                {
                    // loop over all children
                    for ( uint k = 0; k < C; ++k )
                    {
                        // test if child exists
                        if ( mChildren[ k ] != nullptr )
                        {
                            mChildren[ k ]->flag_descendants();
                        }
                    }
                }
            }
        }
        // -----------------------------------------------------------------------------

        void
        unflag_descendants()
        {
            if ( mFlag == true )
            {
                // unflag myself
                mFlag = false;

                // flag children
                // test if children exist
                if ( mChildrenFlag )
                {
                    // loop over all children
                    for ( uint k = 0; k < C; ++k )
                    {
                        // test if child exists
                        if ( mChildren[ k ] != nullptr )
                        {
                            mChildren[ k ]->unflag_descendants();
                        }
                    }
                }
            }
        }

        // -----------------------------------------------------------------------------

        // counts fagged basis
        luint
        count_descendants() override
        {
            // Initialize counter
            luint tBasisCount = 0;

            // test if self has been flagged
            if ( mFlag )
            {
                // add self to counter
                tBasisCount++;

                // test if children exist
                if ( mChildrenFlag )
                {
                    // count children
                    for ( uint iChild = 0; iChild < C; ++iChild )
                    {
                        // test if child exists
                        if ( mChildren[ iChild ] != nullptr )
                        {
                            tBasisCount += mChildren[ iChild ]->count_descendants();
                        }
                    }
                }

                // flag this basis
                mFlag = false;
            }

            return tBasisCount;
        }

        // -----------------------------------------------------------------------------

        // counts inflagged basis
        void
        collect_descendants( Cell< Basis* >& aBasisList,
                luint&                       aBasisCount ) override
        {
            // test if self has been flagged
            if ( !mFlag )
            {
                // add self to list
                aBasisList( aBasisCount++ ) = this;

                // test if children exist
                if ( mChildrenFlag )
                {
                    // add children to list
                    for ( uint k = 0; k < C; ++k )
                    {
                        // test if child exists
                        if ( mChildren[ k ] != nullptr )
                        {
                            mChildren[ k ]->collect_descendants( aBasisList, aBasisCount );
                        }
                    }
                }

                // flag this basis
                mFlag = true;
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Increments the parent counter. Needed for parent identification.
         */
        void
        increment_parent_counter()
        {
            ++mNumberOfParents;
        }
        // -----------------------------------------------------------------------------

        void
        insert_parent( Basis* aParent )
        {
            if ( mParentCounter == 0 )
            {
                // reserve array
                mParents = new Basis*[ mNumberOfParents ];

                // reset array
                for ( uint k = 0; k < mNumberOfParents; ++k )
                {
                    mParents[ k ] = nullptr;
                }
            }
            mParents[ mParentCounter++ ] = aParent;
        }

        // -----------------------------------------------------------------------------

        Basis*
        get_parent( uint aParentNumber )
        {
            return mParents[ aParentNumber ];
        }

        // -----------------------------------------------------------------------------

        uint
        get_number_of_parents()
        {
            return mNumberOfParents;
        }

        // -----------------------------------------------------------------------------

        void
        set_active_index( luint aIndex )
        {
            mActiveIndex = aIndex;
        }

        //------------------------------------------------------------------------------

        luint
        get_active_index()
        {
            return mActiveIndex;
        }

        //------------------------------------------------------------------------------

        void
        get_basis_local_child_inds( Matrix< DDSMat >& aChildren )
        {
            uint tChildren[ C ];

            // count children
            uint tCount = 0;

            if ( mChildrenFlag )
            {
                for ( uint k = 0; k < C; ++k )
                {
                    if ( mChildren[ k ] != nullptr )
                    {
                        if ( mChildren[ k ]->is_active() || mChildren[ k ]->is_refined() )
                        {
                            tChildren[ tCount++ ] = k;
                        }
                    }
                }
            }

            aChildren.set_size( tCount, 1 );

            for ( uint i = 0; i < tCount; ++i )
            {
                aChildren( i ) = tChildren[ i ];
            }
        }

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------

}    // namespace moris::hmr

#endif /* SRC_HMR_CL_HMR_BSPLINE_HPP_ */

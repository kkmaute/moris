/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Vertex.hpp
 *
 */

#pragma once
#include <limits>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "SDF_Tools.hpp"
#include "cl_MTK_Enums.hpp"

#include "cl_MTK_Vertex.hpp"

namespace moris::sdf
{
    // -----------------------------------------------------------------------------

    class Cell;

    class Facet;

    // -----------------------------------------------------------------------------
    /**
     * The sdf vertex is a wrapper around an MTK vertex.
     * It contains a pointer to the MTK vertex and
     * has the ability to flag nodes
     */
    class Vertex
    {
        //! index
        const moris_index mIndex;

        //! flag telling if vertex is inside
        mtk::Mesh_Region mRegion = mtk::Mesh_Region::UNDEFINED;

            //! flag telling if an SDF has been calculated for this vertex
        bool mHasSDF = false;

        bool mIsCandidate = false;

        bool mFlag = true;

        // current node coords
        Matrix< DDRMat > mNodeCoords;
        Matrix< DDRMat > mOriginalNodeCoords;

        real   mSDF;
        Facet* mClosestFacet = nullptr;

        uint mCellCounter = 0;

        Vector< Cell* > mCells;

        Vector< Vertex* > mNeighbors;

        // -----------------------------------------------------------------------------

      public:
        // -----------------------------------------------------------------------------

        /**
         * constructor
         */
        Vertex( const moris_index       aIndex,
                const Matrix< DDRMat >& aNodeCoords );
        // -----------------------------------------------------------------------------

        /**
         * destructor
         */
        ~Vertex()
        {
            mCells.clear();
            mNeighbors.clear();
        };

        // -----------------------------------------------------------------------------

        const Matrix< F31RMat >&
        get_coords() const
        {
            return mNodeCoords;
        }

        // -----------------------------------------------------------------------------

        void
        set_region( const mtk::Mesh_Region aRegion )
        {
            mRegion = aRegion;
        }

        // -----------------------------------------------------------------------------

        mtk::Mesh_Region get_region() const
        {
            return mRegion;
        }

        // -----------------------------------------------------------------------------

        void
        set_candidate_flag()
        {
            mIsCandidate = true;
        }

        // -----------------------------------------------------------------------------

        void
        unset_candidate_flag()
        {
            mIsCandidate = false;
        }

        // -----------------------------------------------------------------------------

        bool
        is_candidate() const
        {
            return mIsCandidate;
        }

        // -----------------------------------------------------------------------------

        void
        set_sdf_flag()
        {
            mHasSDF = true;
        }

        // -----------------------------------------------------------------------------

        void
        unset_sdf_flag()
        {
            mHasSDF = false;
        }

        // -----------------------------------------------------------------------------

        bool
        has_sdf() const
        {
            return mHasSDF;
        }

        // -----------------------------------------------------------------------------

        moris_index
        get_index() const
        {
            return mIndex;
        }

        // -----------------------------------------------------------------------------

        void
        flag()
        {
            mFlag = true;
        }

        // -----------------------------------------------------------------------------

        void
        unflag()
        {
            mFlag = false;
        }

        // -----------------------------------------------------------------------------

        bool
        is_flagged() const
        {
            return mFlag;
        }

        // -----------------------------------------------------------------------------

        void
        reset()
        {
            mHasSDF       = false;
            mIsCandidate  = false;
            mFlag         = true;
            mSDF          = std::numeric_limits< real >::max();
            mClosestFacet = nullptr;
            mIsInside     = false;
        }

        // -----------------------------------------------------------------------------

        void
        update_udf( Facet& aFacet );

        // -----------------------------------------------------------------------------

        void
        increment_cell_counter()
        {
            ++mCellCounter;
        }

        // -----------------------------------------------------------------------------

        void
        init_cell_container()
        {
            mCells.resize( mCellCounter, nullptr );
            mCellCounter = 0;
        }

        // -----------------------------------------------------------------------------

        void
        insert_cell( Cell* aCell );

        // -----------------------------------------------------------------------------

        uint
        get_number_of_cells() const
        {
            return mCellCounter;
        }

        // -----------------------------------------------------------------------------

        Cell*
        get_cell( const uint aIndex )
        {
            return mCells( aIndex );
        }

        // -----------------------------------------------------------------------------

        void
        init_neighbor_container( const uint aNumberOfNeighbors )
        {
            mNeighbors.resize( aNumberOfNeighbors, nullptr );
        }

        // -----------------------------------------------------------------------------

        void
        insert_neighbor( Vertex* aNeighbor, const uint aNeighborIndex )
        {
            mNeighbors( aNeighborIndex ) = aNeighbor;
        }

        // -----------------------------------------------------------------------------

        uint
        get_number_of_neighbors() const
        {
            return mNeighbors.size();
        }

        // -----------------------------------------------------------------------------

        Vertex*
        get_neighbor( const uint aNeighborIndex )
        {
            return mNeighbors( aNeighborIndex );
        }

        // -----------------------------------------------------------------------------

        Facet*
        get_closest_facet()
        {
            return mClosestFacet;
        }

        // -----------------------------------------------------------------------------

        uint
        sweep();

        // -----------------------------------------------------------------------------

        real
        get_sdf() const
        {
            if  ( mRegion == mtk::Mesh_Region::INSIDE )
            {
                return -mSDF;
            }
            else
            {
                return mSDF;
            }
        }

        // -----------------------------------------------------------------------------

        void
        rotate_coords( const Matrix< F33RMat >& aRotationMatrix );

        // -----------------------------------------------------------------------------

        void
        reset_coords();

        // -----------------------------------------------------------------------------
    };

    //-------------------------------------------------------------------------------
}    // namespace moris::sdf

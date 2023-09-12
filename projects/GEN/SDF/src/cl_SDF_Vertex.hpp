/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Vertex.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_VERTEX_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_VERTEX_HPP_
#include <limits>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Vertex.hpp"

namespace moris
{
    namespace sdf
    {
// -----------------------------------------------------------------------------

        class Cell;

        class Triangle;

// -----------------------------------------------------------------------------
        /**
         * The sdf vertex is a wrapper around an MTK vertex.
         * It contains a pointer to the MTK vertex and
         * has the ability to flag nodes
         */
        class Vertex
        {
            //! index
            const moris_index   mIndex;

            //! flag telling if vertex is inside
            bool                mIsInside = false;

            //! flag telling if an SDF has been calculated for this vertex
            bool                mHasSDF = false;

            bool                mIsCandidate = false;

            bool                mFlag = true;

            // current node coords
            Matrix< F31RMat >   mNodeCoords;
            Matrix< F31RMat >   mOriginalNodeCoords;

            real                mSDF;
            Triangle *          mClosestTriangle = nullptr;

            uint                mCellCounter = 0;

            moris::Cell< Cell * > mCells;

            moris::Cell< Vertex * > mNeighbors;

// -----------------------------------------------------------------------------
        public:
// -----------------------------------------------------------------------------

            /**
             * constructor
             */
            Vertex( const moris_index aIndex,
                    const Matrix< DDRMat > & aNodeCoords );
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

            const Matrix< F31RMat > &
            get_coords() const
            {
                return mNodeCoords;
            }

// -----------------------------------------------------------------------------

            void
            set_inside_flag()
            {
                mIsInside = true;
            }

// -----------------------------------------------------------------------------

            void
            unset_inside_flag()
            {
                mIsInside = false;
            }

// -----------------------------------------------------------------------------

            bool
            is_inside() const
            {
                return mIsInside;
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
                mHasSDF = false;
                mIsCandidate = false;
                mFlag = true;
                mSDF =  std::numeric_limits<real>::max();
                mClosestTriangle = nullptr;
                mIsInside = false;
            }

// -----------------------------------------------------------------------------

            void
            update_udf( Triangle *  aTriangle );

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
            insert_cell( Cell * aCell );

// -----------------------------------------------------------------------------

            uint
            get_number_of_cells() const
            {
                return mCellCounter;
            }

// -----------------------------------------------------------------------------

            Cell *
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
            insert_neighbor( Vertex * aNeighbor, const uint aNeighborIndex )
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

            Vertex *
            get_neighbor( const uint aNeighborIndex )
            {
                return mNeighbors( aNeighborIndex );
            }

// -----------------------------------------------------------------------------

            Triangle *
            get_closest_triangle()
            {
                return mClosestTriangle;
            }

// -----------------------------------------------------------------------------

            uint
            sweep();

// -----------------------------------------------------------------------------

            real
            get_sdf() const
            {
                if( mIsInside )
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
            rotate_coords( const Matrix< F33RMat > & aRotationMatrix );

// -----------------------------------------------------------------------------

            void
            reset_coords();

// -----------------------------------------------------------------------------
        };

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_VERTEX_HPP_ */


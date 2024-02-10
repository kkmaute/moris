/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Cell.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_CELL_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_CELL_HPP_

#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_SDF_Vertex.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        /**
         * a wrapper for the MTK cell
         */
        class Cell
        {
            const moris_index mIndex;
            const moris_id    mID;

            // cell with MTK vertices
            Vector< Vertex * > mVertices;

            // flag telling if element is in volume
            bool mElementIsInVolume = false;

            // flag telling if element is in surface
            bool mElementIsOnSurface= false;

            // general purpose flag
            bool mFlag = false;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            Cell(   const moris_index          aIndex,
                    const moris_id             aID,
                    const Matrix< IndexMat > & aIndices,
                    Vector< Vertex * >  & aAllVertices );

//-------------------------------------------------------------------------------

            uint
            get_number_of_vertices() const
            {
                return mVertices.size();
            }

//-------------------------------------------------------------------------------

            Vertex *
            get_vertex( const uint & aIndex )
            {
                return mVertices( aIndex );
            }

//-------------------------------------------------------------------------------

            Vector< Vertex * > &
            get_vertices()
            {
                return mVertices;
            }

//-------------------------------------------------------------------------------
            void
            set_volume_flag()
            {
                mElementIsInVolume = true;
            }

//-------------------------------------------------------------------------------

            void
            unset_volume_flag()
            {
                mElementIsInVolume = false;
            }

//-------------------------------------------------------------------------------

            bool
            is_in_volume() const
            {
                return mElementIsInVolume;
            }

//-------------------------------------------------------------------------------

            void
            set_surface_flag()
            {
                mElementIsOnSurface = true;
            }

//-------------------------------------------------------------------------------

            void
            unset_surface_flag()
            {
                mElementIsOnSurface = false;
            }

//-------------------------------------------------------------------------------

            bool
            is_on_surface() const
            {
                return mElementIsOnSurface;
            }

//-------------------------------------------------------------------------------

            void
            flag()
            {
                mFlag = true;
            }

//-------------------------------------------------------------------------------

            void
            unflag()
            {
                mFlag = false;
            }

//-------------------------------------------------------------------------------

            bool
            is_flagged() const
            {
                return mFlag;
            }

//-------------------------------------------------------------------------------

            moris_index
            get_index()
            {
                return mIndex;
            }

//-------------------------------------------------------------------------------

            moris_id
            get_id()
            {
                return mID;
            }

//-------------------------------------------------------------------------------

            real
            get_buffer_diagonal();

//-------------------------------------------------------------------------------
        };
//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_CELL_HPP_ */


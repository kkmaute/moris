/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Data.hpp
 *
 */

#pragma once

#include "typedefs.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Object.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        struct Data
        {

            //! object information
            uint                          mDimension;
            moris::Cell< Facet* >&        mFacets;
            moris::Cell< Facet_Vertex* >& mVertices;
            const uint                    mNumberOfFacets;    // !< number of facets in object

            //!< counter for unsure nodes in voxelizing algorithm
            uint mUnsureNodesCount;

            //!< matrices useful for intersection computation
            // #ifdef MORIS_USE_ARMA
            //             arma::Mat< real > mFacetMinCoordsX;    //!< min coordinate x of triangle bounding box
            //             arma::Mat< real > mFacetMinCoordsY;    //!< min coordinate y of triangle bounding box
            //             arma::Mat< real > mFacetMinCoordsZ;    //!< min coordinate z of triangle bounding box
            //             arma::Mat< real > mFacetMaxCoordsX;    //!< max coordinate x of triangle bounding box
            //             arma::Mat< real > mFacetMaxCoordsY;    //!< max coordinate y of triangle bounding box
            //             arma::Mat< real > mFacetMaxCoordsZ;    //!< max coordinate x of triangle bounding box

            //             arma::uvec        mCandI;              //!< temporary variable needed for triangle preselection
            //             arma::uvec        mCandJ;              //!< temporary variable needed for triangle preselection
            //             arma::uvec        mCandK;
            // #else


            real mBufferDiagonal;

            // counter for volume elements
            uint mVolumeElements = 0;

            // counter for surface elements
            uint mSurfaceElements = 0;

            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            Data( Object& aObject );

            //-------------------------------------------------------------------------------

            ~Data(){};

            //-------------------------------------------------------------------------------
        };
        //-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

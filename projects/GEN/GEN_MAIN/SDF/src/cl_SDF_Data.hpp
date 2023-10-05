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
            //! cell with triangles
            moris::Cell< Facet* >& mFacets;

            moris::Cell< Facet_Vertex* >& mVertices;

            const uint mNumberOfFacets;    // !< number of facets in object

            //!< counter for unsure nodes in voxelizing algorithm
            uint mUnsureNodesCount;
            // #ifdef MORIS_USE_ARMA
            //             arma::Mat< real > mFacetMinCoordsX;       //!< min coordinate x of triangle bounding box
            //             arma::Mat< real > mFacetMinCoordsY;       //!< min coordinate y of triangle bounding box
            //             arma::Mat< real > mTriangleMinCoordsZ;    //!< min coordinate z of triangle bounding box
            //             arma::Mat< real > mFacetMaxCoordsX;       //!< max coordinate x of triangle bounding box
            //             arma::Mat< real > mFacetMaxCoordsY;       //!< max coordinate y of triangle bounding box
            //             arma::Mat< real > mTriangleMaxCoordsZ;    //!< max coordinate x of triangle bounding box
            //             arma::uvec        mCandI;                 //!< temporary variable needed for triangle preselection
            //             arma::uvec        mCandJ;                 //!< temporary variable needed for triangle preselection
            //             arma::uvec        mCandK;
            // #else
            Matrix< DDRMat > mFacetMinCoords;    //!< min coordinates of the facet bouding box
            Matrix< DDRMat > mFacetMaxCoords;    //!< max coordinates of the facet bouding box
            Matrix< DDUMat > mCandJ;             //!< temporary variable needed for triangle preselection
                                                 // #endif

            Matrix< DDRMat > mCoordsK;    //!< temporary variable needed for voxelizing, coordinates of intersection in triangles by a ray
            Matrix< DDUMat > mCandidateFacets;

            moris::Cell< Facet* > mIntersectedTriangles;

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

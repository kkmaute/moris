/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Data.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_DATA_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_DATA_HPP_

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
            moris::Cell< Facet* >        & mFacets;

            moris::Cell< Facet_Vertex* > & mVertices;

            const uint           mNumberOfTriangles;            // !< number of triangles in object

            //!< counter for unsure nodes in voxelizing algorithm
            uint                 mUnsureNodesCount;
#ifdef MORIS_USE_ARMA
            arma::Mat<real> mTriangleMinCoordsX;    //!< min coordinate x of triangle bounding box
            arma::Mat<real> mTriangleMinCoordsY;    //!< min coordinate y of triangle bounding box
            arma::Mat<real> mTriangleMinCoordsZ;    //!< min coordinate z of triangle bounding box
            arma::Mat<real> mTriangleMaxCoordsX;    //!< max coordinate x of triangle bounding box
            arma::Mat<real> mTriangleMaxCoordsY;    //!< max coordinate y of triangle bounding box
            arma::Mat<real> mTriangleMaxCoordsZ;    //!< max coordinate x of triangle bounding box
            arma::uvec mCandI;                             //!< temporary variable needed for triangle preselection
            arma::uvec mCandJ;                             //!< temporary variable needed for triangle preselection
            arma::uvec mCandK;
#else
            Matrix< DDRMat > mTriangleMinCoordsX;    //!< min coordinate x of triangle bounding box
            Matrix< DDRMat > mTriangleMinCoordsY;    //!< min coordinate y of triangle bounding box
            Matrix< DDRMat > mTriangleMinCoordsZ;    //!< min coordinate z of triangle bounding box
            Matrix< DDRMat > mTriangleMaxCoordsX;    //!< max coordinate x of triangle bounding box
            Matrix< DDRMat > mTriangleMaxCoordsY;    //!< max coordinate y of triangle bounding box
            Matrix< DDRMat > mTriangleMaxCoordsZ;    //!< max coordinate x of triangle bounding box
            Matrix< DDUMat > mCandJ;                 //!< temporary variable needed for triangle preselection
#endif

            Matrix< DDRMat > mCoordsK;                //!< temporary variable needed for voxelizing, coordinates of intersection in triangles by a ray
            Matrix< DDUMat > mCandidateFacets;

            moris::Cell< Facet * > mIntersectedTriangles;

            real mBufferDiagonal;

            // counter for volume elements
            uint mVolumeElements = 0;

            // counter for surface elements
            uint mSurfaceElements = 0;

//-------------------------------------------------------------------------------
        public :
//-------------------------------------------------------------------------------

            Data( Object & aObject );

//-------------------------------------------------------------------------------

            ~Data(){};

//-------------------------------------------------------------------------------
        private:
//-------------------------------------------------------------------------------

            void
            init_triangles ();

        };
//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_DATA_HPP_ */


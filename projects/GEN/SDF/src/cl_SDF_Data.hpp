/*
 * cl_SDF_Data.hpp
 *
 *  Created on: Sep 30, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_DATA_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_DATA_HPP_


#include "typedefs.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Triangle_Mesh.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        struct Data
        {
            //! cell with triangles
            Cell< Triangle * > & mTriangles;

            const uint           mNumberOfTriangles;            // !< number of triangles in object

            //! correction list for voxelizing algorithm
            Matrix< IndexMat >   mUnsureNodes;

            Matrix< IndexMat >   mUnsureNodesNew;

            //!< counter for unsure nodes in voxelizing algorithm
            uint                 mUnsureNewNodesCount;

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

            Matrix< DDUMat > mCandidateTriangles;     //!< triangle candidates to be checked for intersection
            Matrix< DDUMat > mIntersectedTriangles;

            Matrix< DDRMat > mCoordsK;                //!< temporary variable needed for voxelizing

//-------------------------------------------------------------------------------
        public :
//-------------------------------------------------------------------------------

            Data( Triangle_Mesh & aMesh ) :
                mTriangles( aMesh.get_triangles() ),
                mNumberOfTriangles( mTriangles.size() )
            {

            }

//-------------------------------------------------------------------------------

            ~Data(){};

//-------------------------------------------------------------------------------
        };

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */


#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_DATA_HPP_ */

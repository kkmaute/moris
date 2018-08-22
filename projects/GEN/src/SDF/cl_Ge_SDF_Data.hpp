/*
 * cl_GE_SDF_Data.hpp
 *
 *  Created on: Mar 7, 2018
 *      Author: messe
 */

#ifndef SRC_GEOMENG_CL_GE_SDF_DATA_HPP_
#define SRC_GEOMENG_CL_GE_SDF_DATA_HPP_

#ifdef MORIS_USE_ARMA
#include <armadillo>
#else
#include <Eigen/Core>
#endif

#include "chronos.hpp"
#include "cl_Cell.hpp" // CON/src
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Mat.hpp" // LNA/src
#include "GeUtilities.hpp"
#include "cl_Ge_SDF_Triangle.hpp"

namespace ge {

// =============================================================================

    struct SDF_Data
    {

// -----------------------------------------------------------------------------

        /**
          * @brief   Block for user settings
          *
          */
        struct Settings
        {
            moris::uint mCandidateSearchDepth  = 2;
            moris::real mBufferDiagonalEpsilon = 0.01;
        };

// -----------------------------------------------------------------------------

        moris::Mat< moris::real >& mLocalSDF;            // !< SDF in local node coordinates
        const moris::uint mNumberOfTriangles;            // !< number of triangles in object
#ifdef MORIS_USE_ARMA
        arma::Mat<moris::real> mTriangleMinCoordsX;    //!< min coordinate x of triangle bounding box
        arma::Mat<moris::real> mTriangleMinCoordsY;    //!< min coordinate y of triangle bounding box
        arma::Mat<moris::real> mTriangleMinCoordsZ;    //!< min coordinate z of triangle bounding box
        arma::Mat<moris::real> mTriangleMaxCoordsX;    //!< max coordinate x of triangle bounding box
        arma::Mat<moris::real> mTriangleMaxCoordsY;    //!< max coordinate y of triangle bounding box
        arma::Mat<moris::real> mTriangleMaxCoordsZ;    //!< max coordinate x of triangle bounding box
        arma::uvec mCandI;                             //!< temporary variable needed for triangle preselection
        arma::uvec mCandJ;                             //!< temporary variable needed for triangle preselection
        arma::uvec mCandK;
#else
        moris::Mat<moris::real> mTriangleMinCoordsX;    //!< min coordinate x of triangle bounding box
        moris::Mat<moris::real> mTriangleMinCoordsY;    //!< min coordinate y of triangle bounding box
        moris::Mat<moris::real> mTriangleMinCoordsZ;    //!< min coordinate z of triangle bounding box
        moris::Mat<moris::real> mTriangleMaxCoordsX;    //!< max coordinate x of triangle bounding box
        moris::Mat<moris::real> mTriangleMaxCoordsY;    //!< max coordinate y of triangle bounding box
        moris::Mat<moris::real> mTriangleMaxCoordsZ;    //!< max coordinate x of triangle bounding box
        moris::Mat<moris::uint> mCandJ;                 //!< temporary variable needed for triangle preselection
#endif
        moris::Cell< ge::SDF_Triangle > mTriangles;              //!< Cell the triangles of the object

        moris::Mat< moris::uint > mCandidateTriangles;     //!< triangle candidates to be checked for intersection
        moris::Mat< moris::uint > mIntersectedTriangles;   //!< triangle which we really intersected
        moris::Mat< moris::real > mCoordsK;                //!< temporary variable needed for voxelizing
        moris::Mat< moris::uint > mUnsureNodes;            //!< correction list for voxelizing algorithm
        moris::Mat< moris::uint > mUnsureNodesNew;         //!< correction list for voxelizing algorithm
        moris::uint               mUnsureNewNodesCount;    //!< counter for unsure nodes in voxelizing algorithm

        Settings mSettings;                                //!< Struct containing settings
        moris::BoostBitset mLocalNodeInsideFlags;          //!< Flags if a node is inside or outside the domain
        moris::BoostBitset mLocalNodeCandidateFlags;       //!< flags if a node is a candidate for SDF calculation
        moris::BoostBitset mLocalNodeCandidateFlagsOld;    //!< needed if mSettings.mCandidateSearchDepth > 0
        moris::Mat< moris::uint > mLocalElementsAtSurface; //!< Elements which intersect with surface
        moris::Mat< moris::uint > mLocalElementsInVolume;  //!< Elements which are fully inside the volume
        moris::real mBufferDiagonal;                       //!< buffer diagonal for element search
        moris::Mat< moris::uint > mLocalCandidateNodes;    //!< the bitset mLocalNodeInsideFlags converted to a vector
        moris::BoostBitset& mLocalNodeHasSdfFlag;               //!< flags if an SDF was calculated for this node

// -----------------------------------------------------------------------------
    public:
// -----------------------------------------------------------------------------
       SDF_Data(moris::Mat< moris::real >& aLocalSDF,
                moris::BoostBitset& aSDFBitset,
                const moris::Mat<moris::uint> &aTriangleTopology,
                const moris::Mat<moris::real> &aTriangleNodeCoords);

// -----------------------------------------------------------------------------

        /**
        * @brief   Initialize data fields for calculation
        *
        * @param[in] aNumberOfNodes        number of local nodes on proc
        * @param[in] aNumberOfElements     number of local elements on proc
        */
        void
        init_data_fields (
                const moris::uint aNumberOfNodes,
                const moris::uint aNumberOfElements );

// -----------------------------------------------------------------------------

        uint get_number_of_inside_nodes()
        {
            return mLocalNodeCandidateFlags.count();
        }

// -----------------------------------------------------------------------------
    private:
// -----------------------------------------------------------------------------

        /**
         * @brief  subroutine for constructor
         *
         */
        void
        init_triangles (
                const moris::Mat<moris::uint> &aTriangleTopology,
                const moris::Mat<moris::real> &aTriangleNodeCoords);

    };
// =============================================================================
} /* namespace ge */

#endif /* SRC_GEOMENG_CL_GE_SDF_DATA_HPP_ */

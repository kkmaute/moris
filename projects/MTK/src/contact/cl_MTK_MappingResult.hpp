/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_MappingResult.hpp
 *
 */
#ifndef MORIS_CL_MTK_MAPPINGRESULT_HPP
#define MORIS_CL_MTK_MAPPINGRESULT_HPP

#include <ostream>
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"
#include "cl_Json_Object.hpp"


namespace moris::mtk
{
    struct MappingResult
    {

        MappingResult( moris_index aSourceMeshIndex, uint aPhysicalDimension, uint aNumberOfPoints );

        /**
         * @brief The index of the source mesh.
         */
        moris_index mSourceMeshIndex;

        /**
         * @brief The physical coordinates of the point from which the mapping was performed.
         */
        Matrix< DDRMat > mSourcePhysicalCoordinate;

        /**
         * @brief The n-th entry is the cell index of the n-th matching point in the mParametricCoordinates matrix.
         */
        Vector< moris_index > mSourceCellIndex;

        /**
         * @brief The cluster index (w.r.t the side set) of the point from which the mapping was performed.
         */
        Vector< moris_index > mSourceClusterIndex;

        /**
         * @brief The cluster index of the mapped point w.r.t the side set given in mTargetSideSetIndices.
         */
        Vector< moris_index > mTargetClusterIndex;

        /**
         * @brief The n-th entry is the cell index of the n-th matching point in the mParametricCoordinates matrix.
         * If the mapping was not successful, the cell index is -1.
         */
        Vector< moris_index > mTargetCellIndices;

        /**
         * @brief The physical coordinates of the point on which the mapping was performed.
         */
        Matrix< DDRMat > mTargetPhysicalCoordinate;

        /**
         * @brief The parametric coordinates of the point on which the mapping was performed.
         */
        Matrix< DDRMat > mTargetParametricCoordinate;

        /**
         * @brief The n-th entry is the index of the Side_Set index of the n-th matching point in the mParametricCoordinates matrix.
         * If the mapping was not successful, the Side_Set index is -1.
         */
        Vector< moris_index > mTargetSideSetIndices;

        /**
         * @brief Contains the normals that were used for the mapping.
         */
        Matrix< DDRMat > mNormals;

        /**
         * @brief Contains the reference normals (i.e. the normals of the undeformed mesh) that were used for the mapping.
         */
        Matrix<DDRMat > mReferenceNormals;

        /**
         * @brief The n-th entry is the distance of the n-th to the mapped point in physical coordinates.
         */
        Vector< real > mSignedDistance;


        Json to_json();
    };
}    // namespace moris::mtk


#endif    // MORIS_CL_MTK_MAPPINGRESULT_HPP

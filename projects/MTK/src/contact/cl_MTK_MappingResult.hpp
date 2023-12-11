//
// Created by frank on 12/8/23.
//

#ifndef MORIS_CL_MTK_MAPPINGRESULT_HPP
#define MORIS_CL_MTK_MAPPINGRESULT_HPP

#include <ostream>
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "typedefs.hpp"
#include "cl_Json_Object.hpp"


namespace moris::mtk
{
    struct MappingResult
    {

        MappingResult( uint aPhysicalDimension, uint aNumberOfPoints );

        /**
         * @brief The physical coordinates of the point from which the mapping was performed.
         */
        Matrix< DDRMat > mSourcePhysicalCoordinate;

        /**
         * @brief The cluster index (w.r.t the side set) of the point from which the mapping was performed.
         */
        moris::Cell< moris_index > mSourceClusterIndex;

        /**
         * @brief The cluster index of the mapped point w.r.t the side set given in mTargetSideSetIndices.
         */
        moris::Cell<moris_index > mTargetClusterIndex;

        /**
         * @brief The n-th entry is the cell index of the n-th matching point in the mParametricCoordinates matrix.
         * If the mapping was not successful, the cell index is -1.
         */
        moris::Cell< moris_index > mTargetCellIndices;

        /**
         * @brief The physical coordinates of the point on which the mapping was performed.
         */
        Matrix< DDRMat > mTargetPhysicalCoordinate;

        /**
         * @brief The parametric coordinates of the point on which the mapping was performed.
         * @attention The coordinate is only valid if the mapping was successful and a corresponding cell index is not -1!
         */
        Matrix< DDRMat > mTargetParametricCoordinate;

        /**
         * @brief The n-th entry is the index of the Side_Set index of the n-th matching point in the mParametricCoordinates matrix.
         * If the mapping was not successful, the Side_Set index is -1.
         */
        moris::Cell< moris_index > mTargetSideSetIndices;

        /**
         * @brief Contains the normals that were used for the mapping.
         */
        Matrix< DDRMat > mNormal;

        /**
         * @brief The n-th entry is the distance of the n-th to the mapped point in physical coordinates.
         */
        moris::Cell< real > mDistances;


        Json to_json();
    };
}    // namespace moris::mtk


#endif    // MORIS_CL_MTK_MAPPINGRESULT_HPP

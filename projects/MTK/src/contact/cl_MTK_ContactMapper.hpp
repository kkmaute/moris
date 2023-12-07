//
// Created by frank on 11/27/23.
//

#ifndef MORIS_CL_MTK_CONTACTMAPPER_HPP
#define MORIS_CL_MTK_CONTACTMAPPER_HPP

#include <ostream>
#include "cl_MTK_Side_Set.hpp"
#include "cl_Cell.hpp"

namespace moris::mtk
{
    class Integration_Mesh;
    class ContactMapper
    {
      public:
        struct MappingResult
        {
            friend std::ostream &operator<<( std::ostream &aOs, MappingResult const &aPoint );

            /**
             * @brief The parametric coordinates of the point on which the mapping was performed.
             * The coordinate is only valid if the mapping was successful and a corresponding cell index is not -1.
             */
            moris::Matrix< DDRMat > mParametricCoordinates;

            /**
             * @brief The n-th entry of this vector is the row index at which the coordinates of the n-th Side_Set
             * are stored in the mParametricCoordinates matrix.
             */
            moris::Cell< moris_index > mSideSetOffsets;

            /**
             * @brief The n-th entry is the cell index of the n-th matching point in the mParametricCoordinates matrix.
             * If the mapping was not successful, the cell index is -1.
             */
            moris::Cell< moris_index > mCellIndices;

            /**
             * @brief The n-th entry is the index of the Side_Set index of the n-th matching point in the mParametricCoordinates matrix.
             * If the mapping was not successful, the Side_Set index is -1.
             */
            moris::Cell< moris_index > mSideSetIndices;
        };

        ContactMapper(
                Integration_Mesh                                           *aIGMesh,
                moris::Cell< Side_Set * >                                  &aSideSets,
                moris::Cell< std::pair< moris_index, moris_index > > const &aCandidatePairs );

        virtual MappingResult map( Matrix< DDRMat > const &aParametricCoordinate ) = 0;

      protected:
        Integration_Mesh                                    *mIntegrationMeshes;
        moris::Cell< mtk::Side_Set * >                       mSideSets;
        moris::Cell< std::pair< moris_index, moris_index > > mCandidatePairs;
    };
}    // namespace moris::mtk
#endif    // MORIS_CL_MTK_CONTACTMAPPER_HPP

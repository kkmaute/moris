//
// Created by frank on 12/13/23.
//

#ifndef CL_MTK_INTEGRATIONPOINTPAIRS_HPP
#define CL_MTK_INTEGRATIONPOINTPAIRS_HPP

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"

namespace moris::mtk
{
    /**
     * \brief This struct stores information about one or multiple integration points that got mapped from the follower side to the leader side.
     * To uniquely identify the mapping, the information about the follower and leader cell index is stored.
     * Additional information about
     */
    class IntegrationPointPairs
    {
      public:
        IntegrationPointPairs(
                moris_index const         &aFollowerCellIndex,
                Matrix< DDRMat > const    &aFollowerCoordinates,
                moris_index const         &aLeaderCellIndex,
                Matrix< DDRMat > const    &aLeaderCoordinates,
                Vector< real > const &aIntegrationWeights,
                Matrix< DDRMat > const    &aNormals )
                : mFollowerCellIndex( aFollowerCellIndex )
                , mFollowerCoordinates( aFollowerCoordinates )
                , mLeaderCellIndex( aLeaderCellIndex )
                , mLeaderCoordinates( aLeaderCoordinates )
                , mNormals( aNormals )
                , mIntegrationWeights( aIntegrationWeights ){};

        [[nodiscard]] moris_index                get_follower_cell_index() const { return mFollowerCellIndex; }
        [[nodiscard]] const Matrix< DDRMat >    &get_follower_coordinates() const { return mFollowerCoordinates; }
        [[nodiscard]] moris_index                get_leader_cell_index() const { return mLeaderCellIndex; }
        [[nodiscard]] const Matrix< DDRMat >    &get_leader_coordinates() const { return mLeaderCoordinates; }
        [[nodiscard]] const Matrix< DDRMat >    &get_normals() const { return mNormals; }
        [[nodiscard]] const Vector< real > &get_integration_weights() const { return mIntegrationWeights; }

      private:
        /**
         * \brief The cell index of the cell on the follower cell from which the mapping was performed.
         */
        moris_index mFollowerCellIndex;

        /**
         * \brief The parametric coordinates of the integration points on the follower side.
         * \details A (p x n) matrix, where p is the parametric dimension and n is the number of integration points.
         */
        Matrix< DDRMat > mFollowerCoordinates;

        /**
         * \brief The cell index of the cell on the leader cluster to which the mapping was performed.
         */
        moris_index mLeaderCellIndex;

        /**
         * \brief The parametric coordinates of the integration points on the leader side.
         * \details A (p x n) matrix, where p is the parametric dimension and n is the number of integration points.
         */
        Matrix< DDRMat > mLeaderCoordinates;

        /**
         * \brief The normal that was used to perform the mapping from the follower side to the leader side.
         * \details A (d x n) matrix, where d is the physical dimension and n is the number of integration points.
         */
        Matrix< DDRMat > mNormals;

        /**
         * \brief A list of integration points for each of the points in the follower/leader coordinate matrices.
         */
        Vector< real > mIntegrationWeights;
    };
}    // namespace moris::mtk

#endif    // CL_MTK_INTEGRATIONPOINTPAIRS_HPP

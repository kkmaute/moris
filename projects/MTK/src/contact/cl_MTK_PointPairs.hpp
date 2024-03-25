/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_PointPairs.hpp
 *
 */
#pragma once

#include "cl_Matrix_Arma_Dynamic.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"

namespace moris::mtk
{
    /**
     * \brief This struct stores information about one or multiple points that got mapped from the leader side to the follower side.
     * To uniquely identify the mapping, the information about the leader and follower cell index is stored.
     * Additional information about the points e.g. the integration weights (if it is a quadrature point) and the normal that was used to perform the mapping is stored.
     */
    class MappingPointPairs
    {
      public:
        MappingPointPairs()
                : mLeaderCellIndex( -1 )
                , mFollowerCellIndex( -1 ){};

        MappingPointPairs(
                moris_index const      &aLeaderCellIndex,
                Matrix< DDRMat > const &aLeaderCoordinates,
                moris_index const      &aFollowerCellIndex,
                Matrix< DDRMat > const &aFollowerCoordinates,
                Vector< real > const   &aPointDistances,
                Matrix< DDRMat > const &aNormals,
                Matrix< DDRMat > const &aReferenceNormals )
                : mLeaderCellIndex( aLeaderCellIndex )
                , mLeaderCoordinates( aLeaderCoordinates )
                , mFollowerCellIndex( aFollowerCellIndex )
                , mFollowerCoordinates( aFollowerCoordinates )
                , mPointDistances( aPointDistances )
                , mNormals( aNormals )
                , mReferenceNormals( aReferenceNormals )
        {
        }

        virtual ~MappingPointPairs() = default;

        [[nodiscard]] moris_index             get_leader_cell_index() const { return mLeaderCellIndex; }
        [[nodiscard]] const Matrix< DDRMat > &get_leader_coordinates() const { return mLeaderCoordinates; }
        [[nodiscard]] moris_index             get_follower_cell_index() const { return mFollowerCellIndex; }
        [[nodiscard]] const Matrix< DDRMat > &get_follower_coordinates() const { return mFollowerCoordinates; }
        [[nodiscard]] const Matrix< DDRMat > &get_normals() const { return mNormals; }
        [[nodiscard]] const Matrix< DDRMat > &get_reference_normals() const { return mReferenceNormals; }
        [[nodiscard]] const Vector< real >   &get_point_distances() const { return mPointDistances; }

      private:
        /**
         * \brief The cell index of the cell on the leader cell from which the mapping was performed.
         */
        moris_index mLeaderCellIndex;

        /**
         * \brief The parametric coordinates of the integration points on the leader side.
         * \details A (p x n) matrix, where p is the parametric dimension and n is the number of integration points.
         */
        Matrix< DDRMat > mLeaderCoordinates;

        /**
         * \brief The cell index of the cell on the follower cluster to which the mapping was performed.
         */
        moris_index mFollowerCellIndex;

        /**
         * \brief The parametric coordinates of the integration points on the follower side.
         * \details A (p x n) matrix, where p is the parametric dimension and n is the number of integration points.
         */
        Matrix< DDRMat > mFollowerCoordinates;

        /**
         * \brief The distances between the integration points on the leader side and the integration points on the follower side.
         */
        Vector< real > mPointDistances;

        /**
         * \brief The normal that was used to perform the mapping from the leader side to the follower side.
         * \details A (d x n) matrix, where d is the physical dimension and n is the number of integration points.
         */
        Matrix< DDRMat > mNormals;

        /**
         * \brief The reference normal (in the undeformed state) that was used to perform the mapping from the leader side to the follower side.
         * \details A (d x n) matrix, where d is the physical dimension and n is the number of integration points.
         */
        Matrix< DDRMat > mReferenceNormals;
    };

    class IntegrationPointPairs : public MappingPointPairs
    {
      public:
        IntegrationPointPairs() = default;

        IntegrationPointPairs(
                moris_index const      &aLeaderCellIndex,
                Matrix< DDRMat > const &aLeaderCoordinates,
                moris_index const      &aFollowerCellIndex,
                Matrix< DDRMat > const &aFollowerCoordinates,
                Vector< real > const   &aIntegrationWeights,
                Vector< real > const   &aPointDistances,
                Matrix< DDRMat > const &aNormals,
                Matrix< DDRMat > const &aReferenceNormals )
                : MappingPointPairs( aLeaderCellIndex, aLeaderCoordinates, aFollowerCellIndex, aFollowerCoordinates, aPointDistances, aNormals, aReferenceNormals )
                , mIntegrationWeights( aIntegrationWeights ){};

        ~IntegrationPointPairs() override = default;

        [[nodiscard]] const Vector< real > &get_integration_weights() const { return mIntegrationWeights; }

      private:
        /**
         * \brief A list of integration points for each of the points in the leader/follower coordinate matrices.
         */
        Vector< real > mIntegrationWeights;
    };

    class NodalPointPairs : public MappingPointPairs
    {
      public:
        NodalPointPairs() = default;

        NodalPointPairs(
                moris_index const           &aLeaderCellIndex,
                Matrix< DDRMat > const      &aLeaderCoordinates,
                Vector< moris_index > const &aLeaderNodeIndices,
                moris_index const           &aFollowerCellIndex,
                Matrix< DDRMat > const      &aFollowerCoordinates,
                Vector< real > const        &aPointDistances,
                Matrix< DDRMat > const      &aNormals,
                Matrix< DDRMat > const      &aReferenceNormals )
                : MappingPointPairs( aLeaderCellIndex, aLeaderCoordinates, aFollowerCellIndex, aFollowerCoordinates, aPointDistances, aNormals, aReferenceNormals )
                , mLeaderNodeIndices( aLeaderNodeIndices ){};

        ~NodalPointPairs() override = default;

        [[nodiscard]] const Vector< moris_index > &get_leader_node_indices() const { return mLeaderNodeIndices; }

      private:
        /**
         * \brief The indices of the nodes on the leader side.
         * \details A vector of size n, where n is the number of nodes.
         */
        Vector< moris_index > mLeaderNodeIndices;
    };
}    // namespace moris::mtk

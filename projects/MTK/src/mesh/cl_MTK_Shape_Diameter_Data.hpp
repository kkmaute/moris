/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Shape_Diameter_Data.hpp
 *
 */

#pragma once
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace moris::mtk
{
    typedef Vector< std::pair< uint, real > > Intersection_Vector;    // pair of (facet index, distance)

    /**
     * Helper function for applying a Gaussian function to a value, used for weighting ray directions in the shape diameter calculation. The standard deviation controls the spread of the weights.
     */
    real normal_dist( real aValue, real aStdDev );

    // Helper struct that holds the ray directions and weights for the shape diameter calculation. Each vertex has a cone of rays emanating from it, and this struct stores the directions and weights for those rays.
    struct Ray_Cones
    {
        Vector< real > mDirectionWeights;    // Weight for each ray direction
        Vector< real > mTheta;               // Angle from the normal for each ray direction (in radians)
        Vector< real > mPhi;                 // Angle from the vertical for each ray direction (in radians)

        Vector< Matrix< DDRMat > > mRayDirections;    // Each entry in the vector corresponds to a vertex. Each Matrix contains the ray directions as columns.

        /**
         * Constructor to size the data structures
         */
        Ray_Cones( uint aNumFacets, uint aNumRays, uint aDim )
                : mDirectionWeights( aNumRays )
                , mTheta( aNumRays )
                , mPhi( aDim == 3 ? aNumRays : 0 )
                , mRayDirections( aNumFacets, Matrix< DDRMat >( aDim, aNumRays ) )
        {
        }
    };

    class Shape_Diameter_Data
    {
        // Stores the weights and directions for the rays in the cone for each vertex. These are precomputed and stored to avoid redundant calculations during ray casting.
        Ray_Cones mRayCones;

        // Each entry in the vector corresponds to a vertex.
        // Each vertex has a number of rays( aNumPolarRays* aNumAzimuthRays ) associated with it.
        // Each ray only stores its nearest intersection (facet index, distance)
        // Usage: Input( vertex index, ray index ) -> Output( facet index, distance )
        Vector< Intersection_Vector > mDistances;

      public:
        /**
         * Constructor
         */
        Shape_Diameter_Data(
                real                    aConeAngle,
                uint                    aNumPolarRays,
                uint                    aNumAzimuthRays,
                const Matrix< DDRMat >& aNormals,
                uint                    aNumRays );

        void set_distances( Vector< Intersection_Vector >& aDistances );

        const Ray_Cones& get_ray_cones();

        const Vector< Intersection_Vector >& get_distances();

        const Intersection_Vector& get_distances( uint aConeIndex );

        const std::pair< uint, real >& get_distance( uint aConeIndex, uint aRayIndex );

        const Vector< Matrix< DDRMat > >& get_directions();

        /**
         * Gets the direction of a ray for a given vertex and ray index.
         */
        Matrix< DDRMat > get_ray_direction( uint aConeIndex, uint aRayIndex );

        /**
         * Gets the weight associated with a given ray. The weight for each ray is the same for every shape diameter calculation
         */
        real get_ray_weight( uint iRayIndex );

        /**
         * Gets the rotation matrix needed to rotate the normal to the correct direction for a given ray index.
         * This is used for computing sensitivities, where we need to know how the ray direction changes with respect to changes in the normal direction.
         * Note that the rotation matrix is the same for all cones, so no index is needed
         */
        Matrix< DDRMat > get_ray_rotation_matrix( uint aRayIndex );

        uint get_num_cones();

        /**
         * Gets the number of rays per cone
         */
        uint get_num_rays_per_cone();

        /**
         * @brief For each raycast result, determines the nearest intersection that is not part of the originating vertex's facets
         * Used for shape diameter computation
         *
         * @param aAllDistances Raycast results for ray for each vertex - result of a batch raycast
         * @return Vector< Intersection_Vector > Nearest non-trivial intersection for each ray at each vertex. Size: <number of vertices> x <number of rays per vertex>
         */
        void store_nearest_nontrivial_intersections( const Vector< Vector< Intersection_Vector > >& aAllDistances );

      private:
        /**
         * Determines the direction that all rays should be cast in based on the normal direction and the cone parameters
         */
        Ray_Cones build_ray_cone_angles( real aConeAngle, uint aNumPolarRays, uint aNumAzimuthRays, const Matrix< DDRMat >& aNormals ) const;
    };
}    // namespace moris::mtk

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Shape_Diameter_Data.cpp
 *
 */

#include "cl_MTK_Shape_Diameter_Data.hpp"

#include "SDF_Tools.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"

namespace moris::mtk
{
    real normal_dist( real aValue, real aStdDev )
    {
        return ( 1.0 / ( aStdDev * std::sqrt( 2.0 * M_PI ) ) ) * std::exp( -0.5 * std::pow( aValue / aStdDev, 2.0 ) );
    }

    Shape_Diameter_Data::Shape_Diameter_Data(
            real                    aConeAngle,
            uint                    aNumPolarRays,
            uint                    aNumAzimuthRays,
            const Matrix< DDRMat >& aNormals,
            uint                    aNumRays )
            : mRayCones( this->build_ray_cone_angles( aConeAngle, aNumPolarRays, aNumAzimuthRays, aNormals ) )
            , mDistances( aNormals.n_cols(), Intersection_Vector( aNumRays ) )
    {
        // Check that the number of rays matches
        MORIS_ERROR( mRayCones.mRayDirections( 0 ).n_cols() == aNumRays,
                "Shape_Diameter_Data::Constructor - Number of rays does not match between Ray_Cones and input parameter." );

        // Check that all weights are positive
        for ( const auto& tWeight : mRayCones.mDirectionWeights )
        {
            MORIS_ERROR( tWeight > 0.0, "Shape_Diameter_Data::Constructor - All direction weights must be positive." );
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    void Shape_Diameter_Data::set_distances( Vector< Intersection_Vector >& aIntersections )
    {
        mDistances = aIntersections;
    }

    // --------------------------------------------------------------------------------------------------------------

    const Ray_Cones& Shape_Diameter_Data::get_ray_cones()
    {
        return mRayCones;
    }

    // --------------------------------------------------------------------------------------------------------------

    const Vector< Matrix< DDRMat > >& Shape_Diameter_Data::get_directions()
    {
        return mRayCones.mRayDirections;
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Shape_Diameter_Data::get_ray_direction( uint aConeIndex, uint aRayIndex )
    {
        MORIS_ASSERT( aConeIndex < mRayCones.mRayDirections.size(), "get_ray_direction - Cone index out of bounds." );
        MORIS_ASSERT( aRayIndex < mRayCones.mRayDirections( aConeIndex ).n_cols(), "get_ray_direction - Ray index out of bounds." );
        return mRayCones.mRayDirections( aConeIndex ).get_column( aRayIndex );
    }

    // --------------------------------------------------------------------------------------------------------------

    const Vector< Intersection_Vector >& Shape_Diameter_Data::get_distances()
    {
        return mDistances;
    }

    // --------------------------------------------------------------------------------------------------------------

    const Intersection_Vector& Shape_Diameter_Data::get_distances( uint aConeIndex )
    {
        MORIS_ASSERT( aConeIndex < mDistances.size(), "get_distance - Cone index out of bounds" );
        return mDistances( aConeIndex );
    }


    // --------------------------------------------------------------------------------------------------------------

    const std::pair< uint, real >& Shape_Diameter_Data::get_distance( uint aConeIndex, uint aRayIndex )
    {
        MORIS_ASSERT( aRayIndex < this->get_num_rays_per_cone(), "get_distance - Ray index out of bounds" );
        return mDistances( aConeIndex )( aRayIndex );
    }

    // --------------------------------------------------------------------------------------------------------------

    real Shape_Diameter_Data::get_ray_weight( uint iRayIndex )
    {
        MORIS_ASSERT( iRayIndex < mRayCones.mDirectionWeights.size(), "get_ray_weight - Ray index out of bounds." );
        return mRayCones.mDirectionWeights( iRayIndex );
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Shape_Diameter_Data::get_ray_rotation_matrix( uint aRayIndex )
    {
        return sdf::rotation_matrix( mRayCones.mTheta( aRayIndex ) );
    }

    // --------------------------------------------------------------------------------------------------------------

    uint Shape_Diameter_Data::get_num_cones()
    {
        return mDistances.size();
    }

    // --------------------------------------------------------------------------------------------------------------

    uint Shape_Diameter_Data::get_num_rays_per_cone()
    {
        return mRayCones.mDirectionWeights.size();
    }

    // --------------------------------------------------------------------------------------------------------------

    void Shape_Diameter_Data::store_nearest_nontrivial_intersections( const Vector< Vector< Intersection_Vector > >& aAllDistances )
    {
        MORIS_ASSERT( aAllDistances.size() == this->get_num_cones(), "Input intersections size does not match number of cones for shape diameter." );

        // Loop over cones
        for ( uint iF = 0; iF < aAllDistances.size(); iF++ )
        {
            MORIS_ASSERT( aAllDistances( iF ).size() == aAllDistances( 0 ).size(), "Inconsistent number of rays per cone in the input intersections." );

            // Loop over rays for this facet
            for ( uint iR = 0; iR < aAllDistances( iF ).size(); iR++ )
            {
                MORIS_ASSERT( aAllDistances( iF )( iR ).size() > 0, "No intersections found for cone %d, ray %d. Normal may be incorrect.", iF, iR );

                // Loop over the intersections and find the first intersection that is not with the facet the ray came from
                for ( uint iI = 0; iI < aAllDistances( iF )( iR ).size(); iI++ )
                {
                    if ( aAllDistances( iF )( iR )( iI ).first != iF )    // Check that we didn't get an intersection with the originating facet by mistake
                    {
                        mDistances( iF )( iR ) = aAllDistances( iF )( iR )( iI );
                        break;
                    }

                    MORIS_ASSERT( iI < aAllDistances( iF )( iR ).size() - 1, "No non-trivial intersection found for cone %d, ray %d.", iF, iR );
                }
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Ray_Cones Shape_Diameter_Data::build_ray_cone_angles( real aConeAngle, uint aNumPolarRays, uint aNumAzimuthRays, const Matrix< DDRMat >& aNormals ) const
    {
        bool tEvenWeights = false;    // brendan delete temporary for debugging
        bool tUseGaussian = true;     // brendan delete temporary for debugging
        if ( tEvenWeights )
        {
            std::cout << "FACET SHAPE DIAMETER - USING EVEN WEIGHTS\n";
        }
        else if ( tUseGaussian )
        {
            std::cout << "FACET SHAPE DIAMETER - USING GAUSSIAN WEIGHTS\n";
        }
        else
        {
            std::cout << "FACET SHAPE DIAMETER - USING 1/THETA WEIGHTS\n";
        }

        const uint tNumFacets = aNormals.n_cols();
        const uint tDim       = aNormals.n_rows();

        // Initialize return variable
        Ray_Cones tRayCones( tNumFacets, aNumPolarRays * aNumAzimuthRays, tDim );

        switch ( tDim )
        {
            case 2:
            {
                MORIS_ASSERT( aNumPolarRays % 2 == 0, "build_ray_cone_angles - Number of polar rays must be even in 2D." );

                const real tdAlpha = aConeAngle / static_cast< real >( aNumPolarRays );

                Matrix< DDRMat > tRotation;

                for ( uint iF = 0; iF < tNumFacets; iF++ )
                {
                    auto tFacetNormal = aNormals.get_column( iF );

                    for ( uint iR = 0; iR < aNumPolarRays; ++iR )
                    {
                        real tRayAngle;
                        if ( iR == 0 )
                        {
                            // Special case for the first ray
                            tRayAngle = tdAlpha;    // Set the predefined angle
                        }
                        else if ( iR == ( aNumPolarRays / 2 ) )    // For the 16th ray (index 15)
                        {
                            // Special case for the 16th ray
                            tRayAngle = -1.0 * tdAlpha;    // Set the same or another predefined angle
                        }
                        else if ( iR < ( aNumPolarRays / 2 ) )
                        {
                            // Positive angles for the first half of the rays (excluding ray 0)
                            tRayAngle = ( iR + 1 ) * tdAlpha;
                        }
                        else
                        {
                            // Negative angles for the second half of the rays (excluding ray 15)
                            tRayAngle = -1.0 * ( iR - ( aNumPolarRays / 2 ) + 1 ) * tdAlpha;
                        }

                        tRayCones.mTheta( iR ) = tRayAngle * M_PI / 180;    // Convert to radians
                        // tRayCones.mDirectionWeights( iR ) = 1.0 / std::abs( tRayAngle );
                        tRayCones.mDirectionWeights( iR ) = tEvenWeights ? 1.0 : tUseGaussian ? normal_dist( tRayAngle, 10.0 )
                                                                                              : 1.0 / std::abs( tRayAngle );    // brendan delete temporary for debugging

                        tRotation = sdf::rotation_matrix( tRayCones.mTheta( iR ) );

                        // The surface mesh has outward normals. Thus, we need to invert the cone direction to shoot inward
                        Matrix< DDRMat > tConeRay = -tRotation * tFacetNormal;

                        // Check that the ray is opposite of the facet normal
                        MORIS_ASSERT( dot( tConeRay, tFacetNormal ) < 0, "build_ray_cone_angles - Ray direction for cone is not opposite to facet normal." );

                        // Check that the norm is not too small
                        MORIS_ERROR( norm( tConeRay ) >= 1e-9, "build_ray_cone_angles - Ray direction has norm <1e-9. Should be unit." );

                        MORIS_ASSERT( norm( tConeRay ) - 1.0 < 1e-12, "build_ray_cone_angles - Ray direction for cone is not unit." );

                        tRayCones.mRayDirections( iF ).set_column( iR, tConeRay );
                    }
                }

                break;
            }
            case 3:
            {
                real coneAngle    = M_PI / 48;
                real midConeAngle = coneAngle / 2.0;    // Maximum polar angle of the cone (30 degrees)
                real epsilon      = 1e-2;

                for ( uint iF = 0; iF < tNumFacets; iF++ )
                {
                    // Get the vertex position and vertex normal
                    Matrix< DDRMat > e3 = aNormals.get_column( iF );    // Facet normal is the e3 vector

                    // Compute `e1` as a vector orthogonal to `e3` (any vector not collinear with `e3`)
                    Matrix< DDRMat > e1;
                    e1.set_size( 3, 1, 0.0 );
                    if ( fabs( e3( 0, 0 ) ) > fabs( e3( 1, 0 ) ) )
                    {
                        e1( 0, 0 ) = -1.0 * e3( 2, 0 );
                        e1( 2, 0 ) = e3( 0, 0 );    // Choose an arbitrary orthogonal vector
                    }
                    else
                    {
                        e1( 1, 0 ) = -1.0 * e3( 2, 0 );
                        e1( 2, 0 ) = e3( 1, 0 );    // Another option if the x-component is small
                    }
                    e1 = e1 / norm( e1 );    // Normalize `e1`

                    // Compute `e2` as orthogonal to both `e1` and `e3`
                    Matrix< DDRMat > e2 = cross( e3, e1 );    // e2 = e3 x e1

                    // Loop through polar angles (latitude) confined by the coneAngle
                    for ( uint iPolar = 0; iPolar < aNumPolarRays; ++iPolar )
                    {
                        real phi;
                        if ( iPolar < aNumPolarRays / 2 )
                        {
                            // Negative part of φ
                            phi = -midConeAngle + epsilon + ( midConeAngle / ( aNumPolarRays / 2 - 1 ) ) * iPolar;
                        }
                        else
                        {
                            // Positive part of φ, avoiding zero
                            phi = epsilon + ( ( midConeAngle - epsilon ) / ( aNumPolarRays / 2 - 1 ) ) * ( iPolar - aNumPolarRays / 2 );
                        }
                        // Loop through azimuthal angles (longitude)
                        for ( uint iAzimuth = 0; iAzimuth < aNumAzimuthRays; ++iAzimuth )
                        {
                            // Get the index for this ray
                            uint iRayIndex = iAzimuth + iPolar * aNumAzimuthRays;

                            // Compute azimuthal angle θ (from 0 to 2π) and store the spherical coordinates
                            tRayCones.mTheta( iRayIndex ) = ( 2.0 * M_PI / aNumAzimuthRays ) * iAzimuth;
                            tRayCones.mPhi( iRayIndex )   = phi;

                            // Convert spherical coordinates to Cartesian coordinates
                            real x = std::sin( phi ) * std::cos( tRayCones.mTheta( iRayIndex ) );
                            real y = std::sin( phi ) * std::sin( tRayCones.mTheta( iRayIndex ) );
                            real z = std::cos( phi );

                            // Construct the ray in the e1-e2-e3 local coordinates
                            Matrix< DDRMat > tRayDirection = x * e1 + y * e2 + z * e3;

                            // Add the vertex normal direction (e3)
                            Matrix< DDRMat > tConeRay = tRayDirection + e3;

                            // Ensure the ray is aligned with the vertex normal
                            if ( dot( tConeRay, e3 ) > 0 )
                            {
                                tConeRay = 1.0 * tConeRay;
                            }
                            else
                            {
                                tConeRay = -1.0 * tConeRay;
                            }

                            // Numerical adjustment if the norm is too small
                            if ( norm( tConeRay ) < 1e-9 )
                            {
                                tConeRay += 1e-12;
                            }

                            // Store the ray weight (optionally use a different metric if needed)
                            if ( iF == 0 )
                            {
                                tRayCones.mDirectionWeights( iRayIndex ) = 1.0 / std::acos( dot( tConeRay, e3 ) / ( norm( tConeRay ) * norm( e3 ) ) );
                            }
#if MORIS_HAVE_DEBUG
                            else
                            {
                                // Check that the weight is consistent across all facets
                                MORIS_ASSERT( std::abs( tRayCones.mDirectionWeights( iRayIndex ) - ( 1.0 / std::acos( dot( tConeRay, e3 ) / ( norm( tConeRay ) * norm( e3 ) ) ) ) ) < 1e-12,
                                        "Inconsistent ray weight for facet %d, polar %d, azimuth %d",
                                        iF,
                                        iPolar,
                                        iAzimuth );
                            }
#endif
                        }
                    }
                }

                break;
            }
            default:
                MORIS_ERROR( false, "Only 2D-3D implementation" );
                break;
        }

        return tRayCones;
    }

    // --------------------------------------------------------------------------------------------------------------

}    // namespace moris::mtk

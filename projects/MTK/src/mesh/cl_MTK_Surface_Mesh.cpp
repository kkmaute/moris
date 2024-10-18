/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Surface_Mesh.cpp
 *
 */

#include "cl_MTK_Surface_Mesh.hpp"
#include "fn_MTK_QuadraturePointMapper_Ray_ArborX_Details.hpp"

#include <random>
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"
#include "fn_eye.hpp"
#include "fn_trans.hpp"
#include "op_elemwise_mult.hpp"

namespace moris::mtk::arborx
{
    template< typename T >
    T coordinate_to_arborx_point( Matrix< moris::DDRMat > const & aMatrix )
    {
        // handle row vector
        if ( aMatrix.n_cols() != 1 )
        {
            float tZCoord = aMatrix.n_cols() == 3 ? aMatrix( 0, 2 ) : 0.0;
            return { static_cast< float >( aMatrix( 0, 0 ) ), static_cast< float >( aMatrix( 0, 1 ) ), tZCoord };
        }
        // handle column vector
        else
        {
            float tZCoord = aMatrix.n_rows() == 3 ? aMatrix( 2, 0 ) : 0.0;
            return { static_cast< float >( aMatrix( 0, 0 ) ), static_cast< float >( aMatrix( 1, 0 ) ), tZCoord };
        }
    }
}    // namespace moris::mtk::arborx

namespace moris::mtk
{
    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh::Surface_Mesh(
            Matrix< DDRMat >                aVertexCoordinates,
            Vector< Vector< moris_index > > aFacetConnnectivity,
            real                            aIntersectionTolerance )
            : mVertexCoordinates( aVertexCoordinates )
            , mDisplacements( aVertexCoordinates.n_rows(), aVertexCoordinates.n_cols(), 0.0 )
            , mFacetConnectivity( aFacetConnnectivity )
            , mIntersectionTolerance( aIntersectionTolerance )
    {
        // Initialize distortion vectors/matrices
        this->reset_coordinates();

        // Compute the normals of the facets
        this->initialize_facet_normals();

        // Construct the ArborX BVH
        this->construct_bvh();
    }

    void Surface_Mesh::set_all_displacements( const Matrix< DDRMat >& aDisplacements )
    {
        MORIS_ASSERT( aDisplacements.n_rows() == this->get_spatial_dimension(), "Number of vertices in displacement matrix does not match number of vertices in mesh" );
        MORIS_ASSERT( aDisplacements.n_cols() == this->get_number_of_vertices(), "Number of dimensions in displacement matrix does not match number of dimensions in mesh" );

        mDisplacements = aDisplacements;

        // Update the normal vector for all the facets
        this->initialize_facet_normals();

        // Update the bounding volume hierarchy
        this->construct_bvh();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::set_vertex_displacement( const uint aVertexIndex, const Matrix< DDRMat >& aDisplacement )
    {
        MORIS_ASSERT( aDisplacement.n_cols() == 1 or aDisplacement.n_rows() == 1, "Displacement matrix must be a vector" );
        if ( aDisplacement.n_rows() == 1 )
        {
            mDisplacements.set_column( aVertexIndex, trans( aDisplacement ) );
        }
        else
        {
            mDisplacements.set_column( aVertexIndex, aDisplacement );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::get_all_vertex_coordinates() const
    {
        return mVertexCoordinates + mDisplacements;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::get_vertex_coordinates( uint aVertexIndex ) const
    {
        return mVertexCoordinates.get_column( aVertexIndex ) + mDisplacements.get_column( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::get_original_vertex_coordinates( const uint aVertexIndex ) const
    {
        return mVertexCoordinates.get_column( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< moris_index > Surface_Mesh::get_facets_vertex_indices( const uint aFacetIndex ) const
    {
        return mFacetConnectivity( aFacetIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::get_all_vertex_coordinates_of_facet( uint aFacetIndex ) const
    {
        // Get the facets vertex indices
        Vector< moris_index > tVertices    = mFacetConnectivity( aFacetIndex );
        uint                  tNumVertices = tVertices.size();

        // Initialize return matrix
        Matrix< DDRMat > tFacetCoordinates( this->get_spatial_dimension(), tNumVertices );

        // Fill the return matrix with the vertex coordinates
        for ( uint iVertexIndex = 0; iVertexIndex < tNumVertices; iVertexIndex++ )
        {
            tFacetCoordinates.set_column( iVertexIndex, this->get_vertex_coordinates( tVertices( iVertexIndex ) ) );
        }

        return tFacetCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Surface_Mesh::get_spatial_dimension() const
    {
        return mVertexCoordinates.n_rows();
    }

    // --------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Surface_Mesh::get_all_facet_normals() const
    {
        return mFacetNormals;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat > Surface_Mesh::get_facet_normal( const uint aFacetIndex ) const
    {
        return mFacetNormals.get_column( aFacetIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Surface_Mesh::get_number_of_facets() const
    {
        return mFacetConnectivity.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Surface_Mesh::get_number_of_vertices() const
    {
        return mVertexCoordinates.n_cols();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Surface_Mesh::get_intersection_tolerance() const
    {
        return mIntersectionTolerance;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::reset_coordinates()
    {
        // Initialize distortion vectors/matrices
        mDisplacements.set_size( this->get_spatial_dimension(), this->get_number_of_vertices(), 0.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Mesh_Region
    Surface_Mesh::get_region_from_raycast( const Matrix< DDRMat >& aPoint ) const
    {
        // Initialize random ray direction vector
        Matrix< DDRMat > tDirection( this->get_spatial_dimension(), 1 );
        if ( this->get_spatial_dimension() == 2 )
        {
            tDirection = { { 0.6398 }, { -0.4472 } };
        }
        else
        {
            tDirection = { { 0.4990 }, { -0.3534 }, { -0.7912 } };
        }

        // Initialize output
        Mesh_Region tRegion = UNDEFINED;

        real tNorm = norm( tDirection );

        // Cast this ray and get the intersection locations
        Vector< real > tIntersections = this->cast_single_ray_distance_only( aPoint, tDirection );

        // Determine the region based on the number of intersections or if any intersection is close to zero
        tRegion = std::any_of( tIntersections.begin(), tIntersections.end(), [ this, &tNorm ]( real aCoord ) { return std::abs( aCoord ) * tNorm < mIntersectionTolerance; } )
                        ? Mesh_Region::INTERFACE
                        : static_cast< Mesh_Region >( tIntersections.size() % 2 );


        return tRegion;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Mesh_Region >
    Surface_Mesh::batch_get_region_from_raycast( Matrix< DDRMat >& aPoint ) const
    {
        // Getspatial dimension and number of points
        uint tNumPoints = aPoint.n_cols();
        uint tDim       = this->get_spatial_dimension();

        // Initialize random ray direction vector
        Matrix< DDRMat > tDirection( tDim, 1 );
        if ( tDim == 2 )
        {
            tDirection = { { 0.6398 }, { -0.4472 } };
        }
        else
        {
            tDirection = { { 0.4990 }, { -0.3534 }, { -0.7912 } };
        }

        // Initialize output
        Vector< Mesh_Region > tRegions( aPoint.n_cols(), UNDEFINED );

        // Cast all of the rays
        Vector< Vector< Intersection_Vector > > tIntersections = this->cast_batch_of_rays( aPoint, tDirection );

        real tNorm = norm( tDirection );

        // Loop through the intersections and determine the region for each point
        for ( uint iPoint = 0; iPoint < tNumPoints; iPoint++ )
        {
            // Get the intersections for this point FIXME: write a function that returns the distances only to avoid the std::pair overhead
            Intersection_Vector tIntersectionsForPoint = tIntersections( iPoint )( 0 );

            // Store the region for this point
            tRegions( iPoint ) = std::any_of( tIntersectionsForPoint.begin(), tIntersectionsForPoint.end(), [ this, &tNorm ]( std::pair< uint, real > aIntersection ) { return std::abs( aIntersection.second ) * tNorm < mIntersectionTolerance; } )
                                       ? Mesh_Region::INTERFACE
                                       : static_cast< Mesh_Region >( tIntersectionsForPoint.size() % 2 );
        }

        return tRegions;
    }

    //--------------------------------------------------------------------------------------------------------------
    Intersection_Vector
    Surface_Mesh::cast_single_ray(
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection ) const
    {
        // Get the facets that the ray could intersect
        Vector< uint > tCandidateFacets = this->preselect_with_arborx( aPoint, aDirection );

        // Compute the intersection locations for all of the candidates, remove duplicates, and sort
        return this->determine_valid_intersections_from_candidates( aPoint, aDirection, tCandidateFacets );
    }

    //--------------------------------------------------------------------------------------------------------------


    Vector< real >
    Surface_Mesh::cast_single_ray_distance_only(
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection ) const
    {
        // Get the facets that the ray could intersect
        Vector< uint > tCandidateFacets = this->preselect_with_arborx( aPoint, aDirection );

        // Initialize return vector that stores intersections and counter for number of valid intersections
        Vector< real > tIntersections( tCandidateFacets.size() );
        uint           tNumberOfValidIntersections = 0;

        for ( uint iCandidate : tCandidateFacets )
        {
            // Compute the intersection location
            real tIntersection = this->moller_trumbore( iCandidate, aPoint, aDirection );

            // If it is valid, add it to the list
            if ( not std::isnan( tIntersection ) )
            {
                tIntersections( tNumberOfValidIntersections++ ) = tIntersection;
            }
        }

        // Trim output vectors
        tIntersections.resize( tNumberOfValidIntersections );

        if ( tNumberOfValidIntersections != 0 )
        {
            // sort the indices of the array based on the intersection values
            std::sort( tIntersections.begin(), tIntersections.end() );

            // make result unique
            uint tNumberOfUniqueIntersections = 0;

            // initialize return vector
            Vector< real > tUniqueIntersections( tNumberOfValidIntersections );

            // set first entry
            tUniqueIntersections( tNumberOfUniqueIntersections++ ) = tIntersections( 0 );

            // find unique entries
            for ( uint iIntersection = 1; iIntersection < tNumberOfValidIntersections; ++iIntersection )
            {
                if ( std::abs( tIntersections( iIntersection ) - tIntersections( iIntersection - 1 ) ) > mIntersectionTolerance )
                {
                    tUniqueIntersections( tNumberOfUniqueIntersections++ ) = tIntersections( iIntersection );
                }
            }

            // chop vector
            tUniqueIntersections.resize( tNumberOfUniqueIntersections );

            return tUniqueIntersections;
        }
        else
        {
            return tIntersections;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< Intersection_Vector > >
    Surface_Mesh::cast_batch_of_rays(
            Matrix< DDRMat >&           aOrigins,
            Vector< Matrix< DDRMat > >& aDirections ) const
    {
        // Get the number of origins
        uint tNumberOfOrigins = aOrigins.n_cols();

        // Get the facets that the ray could intersect
        Vector< Vector< Vector< uint > > > tCandidateFacets = this->batch_preselect_with_arborx( aOrigins, aDirections );

        // Initialize return variables
        Vector< Vector< Intersection_Vector > > tAllIntersections( tNumberOfOrigins );

        // Loop over every origin
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            // Get the number of directions for this origin
            uint tNumberOfDirections = aDirections( iOrigin ).n_cols();

            // Resize the return vector
            tAllIntersections( iOrigin ).resize( tNumberOfDirections );

            // Loop over every direction (now looping over rays)
            for ( uint iDirection = 0; iDirection < tNumberOfDirections; iDirection++ )
            {
                // Get the candidate facets for this ray
                Vector< uint > tCandidateFacetsForRay = tCandidateFacets( iOrigin )( iDirection );

                // Compute which facets are intersected, remove duplicates, sort, and store in output
                tAllIntersections( iOrigin )( iDirection ) = this->determine_valid_intersections_from_candidates(
                        aOrigins.get_column( iOrigin ),
                        aDirections( iOrigin ).get_column( iDirection ),
                        tCandidateFacetsForRay );
            }
        }

        return tAllIntersections;
    }

    // --------------------------------------------------------------------------------------------------------------

    Vector< Vector< Intersection_Vector > >
    Surface_Mesh::cast_batch_of_rays(
            Matrix< DDRMat >& aOrigins,
            Matrix< DDRMat >& aDirections ) const
    {
        // Get the number of origins
        const uint tNumberOfOrigins    = aOrigins.n_cols();
        const uint tNumberOfDirections = aDirections.n_cols();

        // Get the facets that the ray could intersect
        Vector< Vector< Vector< uint > > > tCandidateFacets = this->batch_preselect_with_arborx( aOrigins, aDirections );

        // Initialize return variables
        Vector< Vector< Intersection_Vector > > tAllIntersections( tNumberOfOrigins );

        // Loop over every origin
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            // Resize the return vector
            tAllIntersections( iOrigin ).resize( tNumberOfDirections );

            // Loop over every direction (now looping over rays)
            for ( uint iDirection = 0; iDirection < tNumberOfDirections; iDirection++ )
            {
                // Get the candidate facets for this ray
                Vector< uint > tCandidateFacetsForRay = tCandidateFacets( iOrigin )( iDirection );

                // Compute which facets are intersected, remove duplicates, sort, and store in output
                tAllIntersections( iOrigin )( iDirection ) = this->determine_valid_intersections_from_candidates(
                        aOrigins.get_column( iOrigin ),
                        aDirections.get_column( iDirection ),
                        tCandidateFacetsForRay );
            }
        }

        return tAllIntersections;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh::moller_trumbore(
            uint                    aFacet,
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection ) const
    {
        switch ( this->get_spatial_dimension() )
        {
            case 2:
            {
                return moller_trumbore_2D( aFacet, aPoint, aDirection );
            }
            case 3:
            {
                return moller_trumbore_3D( aFacet, aPoint, aDirection );
            }
            default:
            {
                MORIS_ERROR( false, "Surface Mesh moller trumbore only implemented for 2D (lines) or 3D (triangles) meshes" );
                return std::numeric_limits< real >::quiet_NaN();
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh::moller_trumbore_2D(
            uint                    aFacet,
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection ) const
    {
        // Get the indices of the vertices for the requested facet
        Vector< moris_index > tVertexIndices = mFacetConnectivity( aFacet );

        // Get the vertex coordinates for the requested facet
        Matrix< DDRMat > tFacetCoordinates = this->get_all_vertex_coordinates_of_facet( aFacet );

        // Build the edge vector
        Matrix< DDRMat > tEdge = tFacetCoordinates.get_column( 1 ) - tFacetCoordinates.get_column( 0 );

        // Vector from the origin of the ray to the origin of the first vertex
        Matrix< DDRMat > tRayToVertex;
        if ( aPoint.n_cols() == 1 )
        {
            tRayToVertex = tFacetCoordinates.get_column( 0 ) - aPoint;
        }
        else
        {
            tRayToVertex = tFacetCoordinates.get_column( 0 ) - trans( aPoint );
        }

        // Get the determinate of the edge and the cast direction
        real tDet = aDirection( 0 ) * tEdge( 1 ) - aDirection( 1 ) * tEdge( 0 );

        // If the determinant is close to zero, the ray is parallel or colinear to the line
        if ( std::abs( tDet ) < mIntersectionTolerance )
        {
            // Check for colinearity
            if ( std::abs( tRayToVertex( 0 ) * aDirection( 1 ) - tRayToVertex( 1 ) * aDirection( 0 ) ) < mIntersectionTolerance )
            {
                // Check if the point is in between the vertices
                if ( ( tFacetCoordinates( 0, 0 ) - aPoint( 0 ) ) * ( tFacetCoordinates( 0, 1 ) - aPoint( 0 ) ) < 0.0 )
                {
                    return 0.0;
                }
                // Check if the facet is in the same direction as the ray (FIXME: technically not necessary if preselection works)
                else if ( ( tFacetCoordinates( 0, 0 ) - aPoint( 0 ) ) * aDirection( 0 ) > 0.0 )
                {
                    // The ray will hit a vertex, return the closer one


                    return ( tFacetCoordinates( 0, 0 ) - aPoint( 0 ) ) < ( tFacetCoordinates( 0, 1 ) - aPoint( 0 ) ) ?                                                            //
                                   ( aPoint.n_cols() == 1 ? norm( tFacetCoordinates.get_column( 0 ) - aPoint ) : norm( tFacetCoordinates.get_column( 0 ) - trans( aPoint ) ) )    //
                                                                                                                     : ( aPoint.n_cols() == 1 ? norm( tFacetCoordinates.get_column( 1 ) - aPoint )
                                                                                                                                              : norm( tFacetCoordinates.get_column( 0 ) - trans( aPoint ) ) );
                }
                // The ray and facet are colinear, but the facet is in the opposite direction of the ray, therefore no intersection
                else
                {
                    return std::numeric_limits< real >::quiet_NaN();
                }
            }
            else
            {
                // Lines are parallel, return nan
                return std::numeric_limits< real >::quiet_NaN();
            }
        }    // end if lines are parallel or colinear

        // Compute the inverse determinant
        real tInverseDeterminant = 1.0 / tDet;

        // Solve the 2D system
        real tDistance = ( tRayToVertex( 0 ) * tEdge( 1 ) - tRayToVertex( 1 ) * tEdge( 0 ) ) * tInverseDeterminant;
        real tU        = ( aDirection( 1 ) * tRayToVertex( 0 ) - aDirection( 0 ) * tRayToVertex( 1 ) ) * tInverseDeterminant;

        // Check if the intersection is within the line segment
        if ( tU < -mIntersectionTolerance or tU > 1.0 + mIntersectionTolerance or tDistance < -mIntersectionTolerance )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }

        return tDistance < -mIntersectionTolerance ? ( std::numeric_limits< real >::quiet_NaN() ) : ( std::abs( tDistance ) < mIntersectionTolerance ? 0.0 : tDistance );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh::moller_trumbore_3D(
            uint                    aFacet,
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection ) const
    {
        // Get the vertex coordinates for the requested facet
        Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( aFacet );

        // Build the edge vectors
        Matrix< DDRMat > tEdge1 = tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 );
        Matrix< DDRMat > tEdge2 = tVertexCoordinates.get_column( 2 ) - tVertexCoordinates.get_column( 0 );

        // Compute the determinant of the edges and the cast direction
        Matrix< DDRMat > tP   = cross( aDirection, tEdge2 );
        real             tDet = dot( tEdge1, tP );

        // If the determinant is close to zero, the ray is parallel to the triangle. Return NaN
        if ( tDet > -mIntersectionTolerance and tDet < mIntersectionTolerance )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }

        // compute the inverse of the determinant
        real tInverseDeterminant = 1.0 / tDet;

        // Compute the vector from the origin to the first vertex
        Matrix< DDRMat > tT;
        if ( aPoint.n_cols() == 1 )
        {
            tT = aPoint - tVertexCoordinates.get_column( 0 );
        }
        else
        {
            tT = trans( aPoint ) - tVertexCoordinates.get_column( 0 );
        }

        // Compute the u parameter
        real tU = dot( tT, tP ) * tInverseDeterminant;

        // If the u parameter is < 0.0 or > 1.0, the intersection is outside the triangle
        if ( tU < 0.0 or tU > 1.0 )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }

        // Compute the vector from the origin to the second vertex
        Matrix< DDRMat > tQ = cross( tT, tEdge1 );

        // Compute the v parameter
        real tV = dot( aDirection, tQ ) * tInverseDeterminant;

        // If the v parameter is < 0.0 or > 1.0, the intersection is outside the triangle
        if ( tV < -mIntersectionTolerance or tU + tV > 1.0 + mIntersectionTolerance )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }

        // Compute the distance from the origin to the intersection point
        real tDistance = dot( tEdge2, tQ ) * tInverseDeterminant;

        // Return nan if the distance is negative, 0.0 if it is close to zero, or the distance otherwise
        return tDistance < -mIntersectionTolerance ? ( std::numeric_limits< real >::quiet_NaN() ) : ( std::abs( tDistance ) < mIntersectionTolerance ? 0.0 : tDistance );
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename MemorySpace, typename ExecutionSpace >
    arborx::QueryRays< MemorySpace > Surface_Mesh::build_arborx_ray_batch(
            ExecutionSpace const &      aExecutionSpace,
            Matrix< DDRMat >&           aOrigins,
            Vector< Matrix< DDRMat > >& aDirections )
    {
        // Get the number of origins, directions, and the total number of rays
        const uint tNumOrigins = aOrigins.n_cols();
        uint       tNumRays    = 0;

        // Calculate the total number of rays
        for ( uint iOrigin = 0; iOrigin < tNumOrigins; iOrigin++ )
        {
            tNumRays += aDirections( iOrigin ).n_cols();
        }

        // Initialize the Kokkos View to store the rays
        Kokkos::View< ArborX::Experimental::Ray*, MemorySpace > tRays( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:rays" ), tNumRays );

        {    // Initialize the rays from the mapping result struct (preinitialized points and normals)
            Kokkos::parallel_for(
                    "initialize_ray_origins",
                    Kokkos::RangePolicy< ExecutionSpace >( aExecutionSpace, 0, tNumOrigins ),
                    KOKKOS_LAMBDA( size_t const iOrigin )    //
                    {
                        // Get the origin point for this group of rays
                        const ArborX::Point tOrigin = arborx::coordinate_to_arborx_point< ArborX::Point >( aOrigins.get_column( iOrigin ) );

                        // Get the number of directions for this origin point
                        const uint tNumDirections = aDirections( iOrigin ).n_cols();

                        // Calculate the starting index for the rays
                        uint tPreviousSize = 0;
                        for ( uint i = 0; i < iOrigin; ++i )
                        {
                            tPreviousSize += aDirections( i ).n_cols();
                        }

                        // Build a ray for every direction
                        Kokkos::parallel_for(
                                "initialize_ray_directions",
                                Kokkos::RangePolicy< ExecutionSpace >( aExecutionSpace, 0, tNumDirections ),
                                KOKKOS_LAMBDA( size_t const iDirection )    //
                                {
                                    size_t tIndex = tPreviousSize + iDirection;

                                    tRays( tIndex ) = ArborX::Experimental::Ray{ tOrigin, arborx::coordinate_to_arborx_point< ArborX::Experimental::Vector >( aDirections( iOrigin ).get_column( iDirection ) ) };
                                } );
                    } );
        }

        return arborx::QueryRays< MemorySpace >{ tRays };
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename MemorySpace, typename ExecutionSpace >
    arborx::QueryRays< MemorySpace > Surface_Mesh::build_arborx_ray_batch(
            ExecutionSpace const & aExecutionSpace,
            Matrix< DDRMat >&      aOrigins,
            Matrix< DDRMat >&      aDirections )
    {
        // Get the number of origins, directions, and the total number of rays
        uint const                                              tNumOrigins    = aOrigins.n_cols();
        uint const                                              tNumDirections = aDirections.n_cols();
        uint const                                              tNumRays       = tNumOrigins * tNumDirections;
        Kokkos::View< ArborX::Experimental::Ray*, MemorySpace > tRays( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:rays" ), tNumRays );
        Kokkos::View< moris_index*, MemorySpace >               tCellIndices( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:cell_indices" ), tNumRays );

        // Initialize the rays from the mapping result struct (preinitialized points and normals)
        Kokkos::parallel_for(
                "initialize_ray_origins",
                Kokkos::RangePolicy< ExecutionSpace >( aExecutionSpace, 0, tNumOrigins ),
                KOKKOS_LAMBDA( size_t const iOrigin )    //
                {
                    // Get the origin point for this group of rays
                    ArborX::Point const tOrigin = arborx::coordinate_to_arborx_point< ArborX::Point >( aOrigins.get_column( iOrigin ) );

                    // Build a ray for every direction
                    Kokkos::parallel_for(
                            "initialize_ray_directions",
                            Kokkos::RangePolicy< ExecutionSpace >( aExecutionSpace, 0, tNumDirections ),
                            KOKKOS_LAMBDA( size_t const iDirection )    //
                            {
                                const size_t tIndex = iOrigin * tNumDirections + iDirection;

                                ArborX::Experimental::Ray tRay = ArborX::Experimental::Ray{ tOrigin, arborx::coordinate_to_arborx_point< ArborX::Experimental::Vector >( aDirections.get_column( iDirection ) ) };
                                tRays( tIndex )                = tRay;
                                tCellIndices( tIndex )         = 0;    // FIXME: Create new struct to remove this dummy cell index for the ray
                            } );
                } );

        return arborx::QueryRays< MemorySpace >{ tRays };
    }

    // --------------------------------------------------------------------------------------------------------------

    Vector< uint > Surface_Mesh::preselect_with_arborx( const Matrix< DDRMat >& aPoint, const Matrix< DDRMat >& aDirection ) const
    {
        using ExecutionSpace = Kokkos::DefaultExecutionSpace;
        using MemorySpace    = ExecutionSpace::memory_space;
        ExecutionSpace tExecutionSpace{};

        // Dummy cell index for the ray
        Kokkos::View< moris_index*, MemorySpace > tCellIndices( Kokkos::view_alloc( tExecutionSpace, Kokkos::WithoutInitializing, "view:cell_indices" ), 1 );
        tCellIndices( 0 ) = 0;

        // Build an ArborX ray from the origin and direction
        Kokkos::View< ArborX::Experimental::Ray*, MemorySpace > tRay( Kokkos::view_alloc( tExecutionSpace, Kokkos::WithoutInitializing, "view:rays" ), 1 );
        tRay( 0 ) = ArborX::Experimental::Ray{ arborx::coordinate_to_arborx_point< ArborX::Point >( aPoint ), arborx::coordinate_to_arborx_point< ArborX::Experimental::Vector >( aDirection ) };

        // Build the query rays
        arborx::QueryRays< MemorySpace > tQueryRays{ tRay };

        // Initialize the results and offsets
        Kokkos::View< arborx::QueryResult*, MemorySpace > tResults( "values", 0 );
        Kokkos::View< int*, MemorySpace >                 tOffsets( "offsets", 0 );

        // Query the BVH for the intersected facets
        mBVH.query( tExecutionSpace, tQueryRays, arborx::RayIntersectionCallback< MemorySpace >{ tQueryRays }, tResults, tOffsets );

        // Return the results as a vector of facet indices
        Vector< uint > tFacetIndices( tResults.extent( 0 ) );
        for ( size_t iIntersection = 0; iIntersection < tResults.extent( 0 ); ++iIntersection )
        {
            tFacetIndices( iIntersection ) = tResults( iIntersection ).mBoxIndex;
        }

        return tFacetIndices;
    }

    // --------------------------------------------------------------------------------------------------------------

    Vector< Vector< Vector< uint > > > Surface_Mesh::batch_preselect_with_arborx(
            Matrix< DDRMat >&           aOrigins,
            Vector< Matrix< DDRMat > >& aDirections ) const
    {
        using ExecutionSpace = Kokkos::DefaultExecutionSpace;
        using MemorySpace    = ExecutionSpace::memory_space;
        ExecutionSpace tExecutionSpace{};

        // Get the number of origins
        const uint tNumOrigins = aOrigins.n_cols();

        // Initialize the return vector
        Vector< Vector< Vector< uint > > > tFacetIndices( tNumOrigins );

        // Initialize the results and offsets
        Kokkos::View< arborx::QueryResult*, MemorySpace > tResults( "values", 0 );
        Kokkos::View< int*, MemorySpace >                 tOffsets( "offsets", 0 );

        // Build the struct that holds the Arborx Ray information from the input matrices
        arborx::QueryRays< MemorySpace > tQueryRays = this->build_arborx_ray_batch< MemorySpace >( tExecutionSpace, aOrigins, aDirections );

        // Query the BVH for the intersected facets
        mBVH.query( tExecutionSpace, tQueryRays, arborx::RayIntersectionCallback< MemorySpace >{ tQueryRays }, tResults, tOffsets );

        if ( tResults.extent( 0 ) > 0 )
        {
            // Initialize counter to index the results
            uint tOffsetIndex = 0;

            // For every ray, get the intersected facets associated with each direction
            for ( uint iOrigin = 0; iOrigin < tNumOrigins; iOrigin++ )
            {
                // Get the number of directions for this origin
                const uint tNumDirections = aDirections( iOrigin ).n_cols();

                // Initialize the return vector for this origin and a temporary vector to store intersection results for this ray
                tFacetIndices( iOrigin ).resize( tNumDirections );
                Vector< uint > tFacetIndicesForRay( tOffsets( tOffsetIndex + 1 ) - tOffsets( tOffsetIndex ) );

                // Loop through the directions for this origin point
                for ( uint iDirection = 0; iDirection < tNumDirections; iDirection++ )
                {
                    // Loop through the results that correspond to this ray
                    for ( moris_index iIntersection = tOffsets( tOffsetIndex ); iIntersection < tOffsets( tOffsetIndex + 1 ); iIntersection++ )
                    {
                        tFacetIndicesForRay( iIntersection - tOffsets( tOffsetIndex ) ) = tResults( iIntersection ).mBoxIndex;
                    }

                    // Add this result to the return
                    tFacetIndices( iOrigin )( iDirection ) = tFacetIndicesForRay;

                    // Increment the offset index
                    tOffsetIndex++;
                }
            }
        }

        return tFacetIndices;
    }

    // --------------------------------------------------------------------------------------------------------------

    Vector< Vector< Vector< uint > > > Surface_Mesh::batch_preselect_with_arborx(
            Matrix< DDRMat >& aOrigins,
            Matrix< DDRMat >& aDirections ) const
    {
        using ExecutionSpace = Kokkos::DefaultExecutionSpace;
        using MemorySpace    = ExecutionSpace::memory_space;
        ExecutionSpace tExecutionSpace{};

        // Get the number of origins
        const uint tNumOrigins    = aOrigins.n_cols();
        const uint tNumDirections = aDirections.n_cols();

        // Initialize the return vector
        Vector< Vector< Vector< uint > > > tFacetIndices( tNumOrigins, Vector< Vector< uint > >( tNumDirections ) );

        // Initialize the results and offsets for the ArborX query
        Kokkos::View< arborx::QueryResult*, MemorySpace > tResults( "values", 0 );
        Kokkos::View< int*, MemorySpace >                 tOffsets( "offsets", 0 );

        // Build the struct that holds the Arborx Ray information from the input matrices
        arborx::QueryRays< MemorySpace > tQueryRays = this->build_arborx_ray_batch< MemorySpace >( tExecutionSpace, aOrigins, aDirections );

        // Query the BVH for the intersected facets
        mBVH.query( tExecutionSpace, tQueryRays, arborx::RayIntersectionCallback< MemorySpace >{ tQueryRays }, tResults, tOffsets );

        if ( tResults.extent( 0 ) > 0 )
        {
            // Initialize counter to index the results
            uint tOffsetIndex = 0;

            // For every ray, get the intersected facets associated with each direction
            for ( uint iOrigin = 0; iOrigin < tNumOrigins; iOrigin++ )
            {
                Vector< uint > tFacetIndicesForRay( tOffsets( tOffsetIndex + 1 ) - tOffsets( tOffsetIndex ) );

                // Loop through the directions for this origin point
                for ( uint iDirection = 0; iDirection < tNumDirections; iDirection++ )
                {
                    // Loop through the results that correspond to this ray
                    for ( moris_index iIntersection = tOffsets( tOffsetIndex ); iIntersection < tOffsets( tOffsetIndex + 1 ); iIntersection++ )
                    {
                        tFacetIndicesForRay( iIntersection - tOffsets( tOffsetIndex ) ) = tResults( iIntersection ).mBoxIndex;
                    }

                    // Add this result to the return
                    tFacetIndices( iOrigin )( iDirection ) = tFacetIndicesForRay;

                    // Increment the offset index
                    tOffsetIndex++;
                }
            }
        }

        return tFacetIndices;
    }

    // --------------------------------------------------------------------------------------------------------------

    Intersection_Vector Surface_Mesh::determine_valid_intersections_from_candidates(
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection,
            const Vector< uint >&   aCandidateFacets ) const
    {
        Intersection_Vector tIntersections( aCandidateFacets.size() );
        uint                tNumberOfValidIntersections = 0;

        for ( uint iCandidate : aCandidateFacets )
        {
            // Compute the intersection location
            real tIntersection = this->moller_trumbore( iCandidate, aPoint, aDirection );

            // If it is valid, add it to the list
            if ( not std::isnan( tIntersection ) )
            {
                tIntersections( tNumberOfValidIntersections ).first    = iCandidate;
                tIntersections( tNumberOfValidIntersections++ ).second = tIntersection;
            }
        }

        if ( tNumberOfValidIntersections != 0 )
        {
            // Trim output vectors
            tIntersections.resize( tNumberOfValidIntersections );

            Intersection_Vector tCleanedIntersections = this->postprocess_raycast_results( tIntersections );

            return tCleanedIntersections;
        }
        else
        {
            return {};
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Intersection_Vector Surface_Mesh::postprocess_raycast_results( Intersection_Vector& aIntersections ) const
    {
        uint tNumIntersections = aIntersections.size();

        // sort the indices of the array based on the intersection values
        std::sort( aIntersections.begin(), aIntersections.end(), []( std::pair< uint, real > i, std::pair< uint, real > j ) {
            return i.second < j.second;
        } );

        // make result unique
        uint tNumberOfUniqueIntersections = 0;

        // initialize return vector
        Intersection_Vector tUniqueIntersections( tNumIntersections );

        // set first entry
        tUniqueIntersections( tNumberOfUniqueIntersections++ ) = aIntersections( 0 );

        // find unique entries
        for ( uint iIntersection = 1; iIntersection < tNumIntersections; ++iIntersection )
        {
            if ( std::abs( aIntersections( iIntersection ).second - aIntersections( iIntersection - 1 ).second ) > mIntersectionTolerance )
            {
                tUniqueIntersections( tNumberOfUniqueIntersections ).second  = aIntersections( iIntersection ).second;
                tUniqueIntersections( tNumberOfUniqueIntersections++ ).first = aIntersections( iIntersection ).first;
            }
        }

        // chop vector and return
        tUniqueIntersections.resize( tNumberOfUniqueIntersections );

        return tUniqueIntersections;
    }

    //-------------------------------------------------------------------------------

    void Surface_Mesh::write_to_file( std::string aFilePath ) const
    {
        // Open file for writing
        std::ofstream tFile;
        tFile.open( aFilePath );

        // Write vertices
        for ( uint iVertex = 0; iVertex < this->get_number_of_vertices(); iVertex++ )
        {
            tFile << "v ";
            for ( uint iDimension = 0; iDimension < this->get_spatial_dimension(); iDimension++ )
            {
                tFile << this->get_vertex_coordinates( iVertex )( iDimension ) << " ";
            }
            tFile << std::endl;
        }

        for ( uint iFacet = 0; iFacet < this->get_number_of_facets(); iFacet++ )
        {
            tFile << "f ";
            const Vector< moris_index > tVertexIndices = this->get_facets_vertex_indices( iFacet );
            for ( uint iVertexIndex = 0; iVertexIndex < tVertexIndices.size(); iVertexIndex++ )
            {
                tFile << tVertexIndices( iVertexIndex ) + 1 << " ";
            }
            tFile << std::endl;
        }

        // close file
        tFile.close();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::set_vertex_coordinates( const uint aVertexIndex, const Matrix< DDRMat >& aCoordinates )
    {
        mVertexCoordinates.set_column( aVertexIndex, aCoordinates );
    }

    // --------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::initialize_facet_normals()
    {
        auto const tNumFacets = static_cast< moris::size_t >( mFacetConnectivity.size() );
        uint const tDim       = this->get_spatial_dimension();

        mFacetNormals.resize( tDim, tNumFacets );

        switch ( tDim )
        {
            case 2:
            {
                for ( moris::size_t iFacetIndex = 0; iFacetIndex < tNumFacets; iFacetIndex++ )
                {
                    Matrix< DDRMat > tNormal( 2, 1 );

                    Matrix< DDRMat > tCoords = this->get_all_vertex_coordinates_of_facet( iFacetIndex );

                    // { { tY2 - tY1  }, { tX1 - tX2 } }
                    tNormal( 0 ) = tCoords( 1, 1 ) - tCoords( 1, 0 );
                    tNormal( 1 ) = tCoords( 0, 0 ) - tCoords( 0, 1 );
                    tNormal      = tNormal / norm( tNormal );

                    mFacetNormals.set_column( iFacetIndex, tNormal );
                }
                break;
            }
            case 3:
            {
                for ( moris::size_t iFacetIndex = 0; iFacetIndex < tNumFacets; iFacetIndex++ )
                {
                    Matrix< DDRMat > tCoords = this->get_all_vertex_coordinates_of_facet( iFacetIndex );

                    Matrix< DDRMat > tNormal = cross( tCoords.get_column( 1 ) - tCoords.get_column( 0 ), tCoords.get_column( 2 ) - tCoords.get_column( 0 ) );
                    tNormal                  = tNormal / norm( tNormal );

                    mFacetNormals.set_column( iFacetIndex, tNormal );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Surface Mesh facet normals only implemented for 2D (lines) or 3D (triangles) meshes" );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::construct_bvh()
    {
        // Create dummy pair for surface mesh
        Vector< std::pair< moris_index, Surface_Mesh > > tSurfaceMesh = { std::pair< moris_index, Surface_Mesh >( { 0, *this } ) };

        // Construct the ArborX boxes from this mesh
        ExecutionSpace                    tExecutionSpace{};
        arborx::QueryBoxes< MemorySpace > tQueryBoxes = arborx::construct_query_boxes< MemorySpace >( tExecutionSpace, tSurfaceMesh );

        // Build the bounding volume hierarchy from the boxes
        mBVH = ArborX::BVH< MemorySpace >( tExecutionSpace, tQueryBoxes );
    }
}    // namespace moris::mtk
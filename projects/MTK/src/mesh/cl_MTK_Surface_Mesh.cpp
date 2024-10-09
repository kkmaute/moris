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

#include <random>
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"
#include "fn_eye.hpp"
#include "fn_trans.hpp"
#include "op_elemwise_mult.hpp"


// BRENDAN FORMAT

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
        mDisplacements.set_column( aVertexIndex, aDisplacement );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::append_vertex_displacement( const uint aVertexIndex, const Matrix< DDRMat >& aDisplacement )
    {
        mDisplacements.set_column( aVertexIndex, mDisplacements.get_column( aVertexIndex ) + aDisplacement );
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

    Vector< moris_index > Surface_Mesh::get_facets_vertex_indices( const uint aFacetIndex ) const
    {
        return mFacetConnectivity( aFacetIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::get_all_vertex_coordinates_of_facet( uint aFacetIndex ) const
    {
        // Get the facets vertex indices
        Vector< moris_index > tVertices = mFacetConnectivity( aFacetIndex );

        // Initialize return matrix
        Matrix< DDRMat > tFacetCoordinates( this->get_spatial_dimension(), tVertices.size() );

        // Fill the return matrix with the vertex coordinates
        for ( uint iVertexIndex = 0; iVertexIndex < tVertices.size(); iVertexIndex++ )
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

    Matrix< DDRMat > Surface_Mesh::get_all_facet_normals() const
    {
        return mFacetNormals;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::get_facet_normal( const uint aFacetIndex ) const
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

    void Surface_Mesh::reset_coordinates()
    {
        // Initialize distortion vectors/matrices
        mDisplacements.set_size( this->get_spatial_dimension(), this->get_number_of_vertices(), 0.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Mesh_Region
    Surface_Mesh::get_region_from_raycast( const Matrix< DDRMat >& aPoint ) const
    {
        // Build a random direction vector
        Matrix< DDRMat > tDirection( this->get_spatial_dimension(), 1 );
        switch ( this->get_spatial_dimension() )
        {
            case 2:
            {
                tDirection( 0, 0 ) = 0.7986;
                tDirection( 1, 0 ) = 0.6018;
            }
            case 3:
            {
                tDirection( 0, 0 ) = 0.5642;
                tDirection( 1, 0 ) = 0.4250;
                tDirection( 2, 0 ) = -0.4367;
            }
            default:
            {
                MORIS_ERROR( false, "Surface Mesh raycast only implemented for 2D (lines) or 3D (triangles) meshes" );
                return Mesh_Region::UNDEFINED;
            }
        }

        // Cast this ray and get the intersection locations
        Intersection_Vector tIntersections = this->cast_single_ray( aPoint, tDirection );

        // Determine the region based on the number of intersections or if any intersection is close to zero
        return std::any_of( tIntersections.begin(), tIntersections.end(), [ this ]( std::pair< uint, real > aIntersection ) { return std::abs( aIntersection.second ) < mIntersectionTolerance; } )
                     ? Mesh_Region::INTERFACE
                     : static_cast< Mesh_Region >( tIntersections.size() % 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Mesh_Region >
    Surface_Mesh::batch_get_region_from_raycast( Matrix< DDRMat >& aPoints ) const
    {
        // Get the number of origins
        uint tNumberOfOrigins = aPoints.n_cols();

        // Build a random direction vector
        Matrix< DDRMat > tDirection( this->get_spatial_dimension(), 1 );
        switch ( this->get_spatial_dimension() )
        {
            case 2:
            {
                tDirection( 0, 0 ) = 0.7986;
                tDirection( 1, 0 ) = 0.6018;
            }
            case 3:
            {
                tDirection( 0, 0 ) = 0.5642;
                tDirection( 1, 0 ) = 0.4250;
                tDirection( 2, 0 ) = -0.4367;
            }
            default:
            {
                MORIS_ERROR( false, "Surface Mesh raycast only implemented for 2D (lines) or 3D (triangles) meshes" );
                return Mesh_Region::UNDEFINED;
            }
        }

        // Cast all of the rays
        Vector< Vector< Intersection_Vector > > tIntersections = this->cast_batch_of_rays( aPoints, tDirection );

        // Initialize return vector
        Vector< Mesh_Region > tRegions( aPoints.n_cols() );

        // Loop over all the origins
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            // Determine the region based on the number of intersections or if any intersection is close to zero
            tRegions( iOrigin ) = std::any_of( tIntersections( iOrigin )( 0 ).begin(), tIntersections( iOrigin )( 0 ).end(), [ this ]( std::pair< uint, real > aIntersection ) { return std::abs( aIntersection.second ) < mIntersectionTolerance; } )
                                        ? Mesh_Region::INTERFACE
                                        : static_cast< Mesh_Region >( tIntersections( iOrigin )( 0 ).size() % 2 );
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

        // Initialize return vector that stores intersections and counter for number of valid intersections
        Intersection_Vector tIntersections( tCandidateFacets.size() );
        uint                tNumberOfValidIntersections = 0;

        for ( uint iCandidate : tCandidateFacets )
        {
            // Compute the intersection location
            real tIntersection = this->moller_trumbore( iCandidate, aPoint, aDirection );

            // If it is valid, add it to the list
            if ( not std::isnan( tIntersection ) )
            {
                tIntersections( tNumberOfValidIntersections ).second  = tIntersection;
                tIntersections( tNumberOfValidIntersections++ ).first = iCandidate;
            }
        }

        // Remove duplicates and sort the intersections, return
        return postprocess_raycast_output( tIntersections );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< Intersection_Vector > >
    Surface_Mesh::cast_batch_of_rays(
            Matrix< DDRMat >& aPoints,
            Matrix< DDRMat >& aDirections ) const
    {
        // Get the number of origins and directions
        uint tNumberOfOrigins    = aPoints.n_cols();
        uint tNumberOfDirections = aDirections.n_cols();

        // Initialize return vector
        Vector< Vector< Intersection_Vector > > tIntersections( aPoints.n_cols(), Vector< Intersection_Vector >( aDirections.n_cols() ) );

        // Get all of the candidate facets
        Vector< Vector< Vector< uint > > > tCandidateFacets = this->batch_preselect_with_arborx( aPoints, aDirections );

        // Loop over all the origins
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            for ( uint iDirection = 0; iDirection < tNumberOfDirections; iDirection++ )
            {
                // Initialize the intersection vector for this origin and direction
                Intersection_Vector tIntersectionsForRay( tCandidateFacets( iOrigin )( iDirection ).size() );
                uint                tNumberOfIntersections = 0;

                // Iterate over all the candidate facets
                for ( uint iCandidate : tCandidateFacets( iOrigin )( iDirection ) )
                {
                    // Compute the intersection location
                    real tIntersection = this->moller_trumbore( iCandidate, aPoints.get_column( iOrigin ), aDirections.get_column( iDirection ) );

                    // If it is valid, add it to the list
                    if ( not std::isnan( tIntersection ) )
                    {
                        tIntersectionsForRay( tNumberOfIntersections ).second  = tIntersection;
                        tIntersectionsForRay( tNumberOfIntersections++ ).first = iCandidate;
                    }

                    // Trim the output vector
                    tIntersectionsForRay.resize( tNumberOfIntersections );

                    // Remove duplicates and sort the intersections, store in output vector
                    tIntersections( iOrigin )( iDirection ) = postprocess_raycast_output( tIntersectionsForRay );
                }
            }
        }

        return tIntersections;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< Intersection_Vector > >
    Surface_Mesh::cast_batch_of_rays(
            Matrix< DDRMat >&           aPoints,
            Vector< Matrix< DDRMat > >& aDirections ) const
    {
        // Get the number of origins and directions
        uint tNumberOfOrigins = aPoints.n_cols();

        // Initialize return vector
        Vector< Vector< Intersection_Vector > > tIntersections( aPoints.n_cols() );

        // Get all of the candidate facets
        Vector< Vector< Vector< uint > > > tCandidateFacets = this->batch_preselect_with_arborx( aPoints, aDirections );

        // Loop over all the origins
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            uint tNumberOfDirections = aDirections( iOrigin ).n_cols();
            tIntersections( iOrigin ).resize( tNumberOfDirections );

            for ( uint iDirection = 0; iDirection < tNumberOfDirections; iDirection++ )
            {
                // Initialize the intersection vector for this origin and direction
                Intersection_Vector tIntersectionsForRay( tCandidateFacets( iOrigin )( iDirection ).size() );
                uint                tNumberOfIntersections = 0;

                // Iterate over all the candidate facets
                for ( uint iCandidate : tCandidateFacets( iOrigin )( iDirection ) )
                {
                    // Compute the intersection location
                    real tIntersection = this->moller_trumbore( iCandidate, aPoints.get_column( iOrigin ), aDirections( iOrigin ).get_column( iDirection ) );

                    // If it is valid, add it to the list
                    if ( not std::isnan( tIntersection ) )
                    {
                        tIntersectionsForRay( tNumberOfIntersections ).second  = tIntersection;
                        tIntersectionsForRay( tNumberOfIntersections++ ).first = iCandidate;
                    }

                    // Trim the output vector
                    tIntersectionsForRay.resize( tNumberOfIntersections );

                    // Remove duplicates and sort the intersections, store in output vector
                    tIntersections( iOrigin )( iDirection ) = postprocess_raycast_output( tIntersectionsForRay );
                }
            }
        }

        return tIntersections;
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
        Matrix< DDRMat > tRayToVertex = tFacetCoordinates.get_column( 0 ) - aPoint;

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
                    return ( tFacetCoordinates( 0, 0 ) - aPoint( 0 ) ) < ( tFacetCoordinates( 0, 1 ) - aPoint( 0 ) ) ? norm( tFacetCoordinates.get_column( 0 ) - aPoint ) : norm( tFacetCoordinates.get_column( 1 ) - aPoint );
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
        }

        // Compute the inverse determinant
        real tInverseDeterminant = 1.0 / tDet;

        // Solve the 2D system
        real tDistance = ( tRayToVertex( 0 ) * tEdge( 1 ) - tRayToVertex( 1 ) * tEdge( 0 ) ) * tInverseDeterminant;
        real tU        = ( aDirection( 0 ) * tRayToVertex( 1 ) - aDirection( 1 ) * tRayToVertex( 0 ) ) * tInverseDeterminant;

        // Check if the intersection is within the line segment
        if ( tU < -mIntersectionTolerance or tU > 1.0 + mIntersectionTolerance or tDistance < 0 )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }
        // check if the ray originates from on the facet
        else if ( std::abs( tDistance ) < mIntersectionTolerance )
        {
            // Snap to the facet
            return 0.0;
        }
        else
        {
            return norm( tDistance * aDirection );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Surface_Mesh::moller_trumbore_3D(
            uint                    aFacet,
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection ) const
    {
        // Get the indices of the vertices for the requested facet
        Vector< moris_index > tVertexIndices = mFacetConnectivity( aFacet );

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
            tT = aPoint - tVertexCoordinates.get_column( tVertexIndices( 0 ) );
        }
        else
        {
            tT = trans( aPoint ) - tVertexCoordinates.get_column( tVertexIndices( 0 ) );
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

        // Return the distance
        return norm( tDistance * aDirection );
    }

    //--------------------------------------------------------------------------------------------------------------


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
        arborx::QueryRays< MemorySpace > tQueryRays{ tRay, tCellIndices };

        // Initialize the results and offsets
        Kokkos::View< arborx::QueryResult*, MemorySpace > tResults( "values", 0 );
        Kokkos::View< int*, MemorySpace >                 tOffsets( "offsets", 0 );

        // Query the BVH for the intersected facets
        mBVH.query( tExecutionSpace, tQueryRays, arborx::IntersectionCallback< MemorySpace >{ tQueryRays }, tResults, tOffsets );

        // Return the results as a vector of facet indices
        Vector< uint > tFacetIndices( tResults.extent( 0 ) );
        for ( size_t iIntersection = 0; iIntersection < tResults.extent( 0 ); ++iIntersection )
        {
            tFacetIndices( iIntersection ) = tResults( iIntersection ).mBoxIndex;
        }

        return tFacetIndices;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< Vector< uint > > > Surface_Mesh::batch_preselect_with_arborx(
            Matrix< DDRMat >& aPoints,
            Matrix< DDRMat >& aDirections ) const
    {
        // Initialize Kokkos execution space
        using ExecutionSpace = Kokkos::DefaultExecutionSpace;
        using MemorySpace    = ExecutionSpace::memory_space;
        ExecutionSpace tExecutionSpace{};

        // Get the number of origins and directions
        uint tNumberOfOrigins    = aPoints.n_cols();
        uint tNumberOfDirections = aDirections.n_cols();

        // Initialize the output vector
        Vector< Vector< Vector< uint > > > tFacetIndices( aPoints.n_cols(), Vector< Vector< uint > >( aDirections.n_cols() ) );

        // Construct input struct for ArborX
        arborx::QueryRays< MemorySpace > tQueryRays = arborx::construct_query_rays_from_primitives< MemorySpace >( tExecutionSpace, aPoints, aDirections );

        // Initialize the results and offsets
        Kokkos::View< arborx::QueryResult*, MemorySpace > tResults( "values", 0 );
        Kokkos::View< int*, MemorySpace >                 tOffsets( "offsets", 0 );

        // Query the BVH for the intersected facets
        mBVH.query( tExecutionSpace, tQueryRays, arborx::IntersectionCallback< MemorySpace >{ tQueryRays }, tResults, tOffsets );

        // Loop through the resutls and store the facet indices
        uint tRayIndex = 0;
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            for ( uint iDirection = 0; iDirection < tNumberOfDirections; iDirection++ )
            {
                // Get the number of intersections for this ray
                uint tNumberOfIntersections = tOffsets( tRayIndex + 1 ) - tOffsets( tRayIndex );

                // Initialize the facet indices vector
                Vector< uint > tFacetIndicesForRay( tNumberOfIntersections );

                // Store the facet indices
                for ( uint iIntersection = 0; iIntersection < tNumberOfIntersections; iIntersection++ )
                {
                    tFacetIndicesForRay( iIntersection ) = tResults( tOffsets( tRayIndex ) + iIntersection ).mBoxIndex;
                }

                // Store the facet indices for this ray
                tFacetIndices( iOrigin )( iDirection ) = tFacetIndicesForRay;

                // Increment the ray index
                tRayIndex++;
            }
        }

        return tFacetIndices;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< Vector< uint > > > Surface_Mesh::batch_preselect_with_arborx(
            Matrix< DDRMat >&           aPoints,
            Vector< Matrix< DDRMat > >& aDirections ) const
    {
        // Initialize Kokkos execution space
        using ExecutionSpace = Kokkos::DefaultExecutionSpace;
        using MemorySpace    = ExecutionSpace::memory_space;
        ExecutionSpace tExecutionSpace{};

        // Get the number of origins and directions
        uint tNumberOfOrigins = aPoints.n_cols();

        // Initialize the output vector
        Vector< Vector< Vector< uint > > > tFacetIndices( aPoints.n_cols() );

        // Construct input struct for ArborX
        arborx::QueryRays< MemorySpace > tQueryRays = arborx::construct_query_rays_from_primitives< MemorySpace >( tExecutionSpace, aPoints, aDirections );

        // Initialize the results and offsets
        Kokkos::View< arborx::QueryResult*, MemorySpace > tResults( "values", 0 );
        Kokkos::View< int*, MemorySpace >                 tOffsets( "offsets", 0 );

        // Query the BVH for the intersected facets
        mBVH.query( tExecutionSpace, tQueryRays, arborx::IntersectionCallback< MemorySpace >{ tQueryRays }, tResults, tOffsets );

        // Loop through the resutls and store the facet indices
        uint tRayIndex = 0;
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            uint tNumberOfDirections = aDirections( iOrigin ).n_cols();
            tFacetIndices( iOrigin ).resize( tNumberOfDirections );

            for ( uint iDirection = 0; iDirection < tNumberOfDirections; iDirection++ )
            {
                // Get the number of intersections for this ray
                uint tNumberOfIntersections = tOffsets( tRayIndex + 1 ) - tOffsets( tRayIndex );

                // Initialize the facet indices vector
                Vector< uint > tFacetIndicesForRay( tNumberOfIntersections );

                // Store the facet indices
                for ( uint iIntersection = 0; iIntersection < tNumberOfIntersections; iIntersection++ )
                {
                    tFacetIndicesForRay( iIntersection ) = tResults( tOffsets( tRayIndex ) + iIntersection ).mBoxIndex;
                }

                // Store the facet indices for this ray
                tFacetIndices( iOrigin )( iDirection ) = tFacetIndicesForRay;

                // Increment the ray index
                tRayIndex++;
            }
        }

        return tFacetIndices;
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Vector
    Surface_Mesh::postprocess_raycast_output( Intersection_Vector& aIntersections ) const
    {
        // sort the input based on the distance from the facet
        std::sort( aIntersections.begin(), aIntersections.end(), []( std::pair< uint, real > a, std::pair< uint, real > b ) { return a.second < b.second; } );

        // Initialize the output vector and store the first intersection if needed
        Intersection_Vector tCleanedIntersections( aIntersections.size() );
        uint                tNumberOfUniqueIntersections = 0;
        if ( aIntersections.size() > 0 )
        {
            tCleanedIntersections( tNumberOfUniqueIntersections++ ) = aIntersections( 0 );

            // loop through the rest of the intersections and test if the next intersection is close to the previous one
            for ( uint iIntersection = 1; iIntersection < aIntersections.size(); iIntersection++ )
            {
                // if the next intersection is not close to the previous one, add it to the output
                if ( std::abs( aIntersections( iIntersection ).second - aIntersections( iIntersection - 1 ).second ) > 10.0 * mIntersectionTolerance )
                {
                    tCleanedIntersections( tNumberOfUniqueIntersections++ ) = aIntersections( iIntersection );
                }
            }
        }

        // Trim the output vector
        tCleanedIntersections.resize( tNumberOfUniqueIntersections );

        return tCleanedIntersections;
    }

    //--------------------------------------------------------------------------------------------------------------


    void Surface_Mesh::write_to_file( std::string aFilePath )
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
                tFile << mVertexCoordinates( iDimension, iVertex ) << " ";
            }
            tFile << std::endl;
        }

        // Write facets
        for ( uint iFacet = 0; iFacet < this->get_number_of_facets(); iFacet++ )
        {
            tFile << "f ";
            Vector< moris_index > tIndices = this->get_facets_vertex_indices( iFacet );
            for ( uint iDimension = 0; iDimension < this->get_spatial_dimension(); iDimension++ )
            {
                tFile << tIndices( iDimension ) + 1 << " ";
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

                    Matrix< DDRMat > tNormal = cross( tCoords.get_column( 2 ) - tCoords.get_column( 0 ), tCoords.get_column( 1 ) - tCoords.get_column( 0 ) );

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

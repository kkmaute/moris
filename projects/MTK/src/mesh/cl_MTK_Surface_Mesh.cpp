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
#include "fn_trans.hpp"
#include "op_elemwise_mult.hpp"
#include <set>

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
    // helper function for 2d raycast
    real cross_2d( const Matrix< DDRMat >& aVector1, const Matrix< DDRMat >& aVector2 )
    {
        return aVector1( 0 ) * aVector2( 1 ) - aVector1( 1 ) * aVector2( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh::Surface_Mesh(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnnectivity,
            real                                   aIntersectionTolerance )
            : mVertexCoordinates( aVertexCoordinates )
            , mDisplacements( aVertexCoordinates.n_rows(), aVertexCoordinates.n_cols(), 0.0 )
            , mFacetConnectivity( aFacetConnnectivity )
            , mIntersectionTolerance( aIntersectionTolerance )
    {
        // Remove the unused vertices, update the facet connectivity
        this->clean_extraneous_vertices();

        // Initialize distortion vectors/matrices
        this->reset_coordinates();

        // Compute the normals of the facets
        this->initialize_facet_normals();

        // Build the vertex connectivity
        this->build_vertex_connectivity();

#if MORIS_HAVE_ARBORX
        // Construct the ArborX BVH
        this->construct_bvh();
#else
        MORIS_LOG_WARNING( "You are using an mtk::Surface_Mesh without ArborX turned on. While all functionality is available, raycasting will be MUCH slower than you'd like." );
#endif
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::clean_extraneous_vertices()
    {
        // STEP 1: Check if all vertices are used in the mesh
        // Initialize vector to store the indices of all the vertices that are connected to facets
        Vector< moris_index > tUsedVertexIndices( mFacetConnectivity.size() * mFacetConnectivity( 0 ).size(), MORIS_INDEX_MAX );

        // Loop through all the facets and store the vertices they use
        uint tIndex = 0;
        for ( auto tFacet : mFacetConnectivity )
        {
            for ( moris_index tVertexIndex : tFacet )
            {
                tUsedVertexIndices( tIndex++ ) = tVertexIndex;
            }
        }

        // Remove duplicates, sort in ascending order
        std::sort( tUsedVertexIndices.begin(), tUsedVertexIndices.end() );
        auto tLast = std::unique( tUsedVertexIndices.begin(), tUsedVertexIndices.end() );
        tUsedVertexIndices.resize( std::distance( tUsedVertexIndices.begin(), tLast ) );

        // STEP 2: Update the vertex coordinates and facet connectivity
        if ( tUsedVertexIndices.size() < mVertexCoordinates.n_cols() )
        {
            // Create map to track the original vertex index to the cleaned vertex index
            std::map< moris_index, moris_index > tVertexIndexMap;

            // Get the number of used vertices
            uint tNumUsedVertices = tUsedVertexIndices.size();

            // Loop through the used vertices update the vertex coordinates
            for ( uint iVertex = 0; iVertex < tNumUsedVertices; ++iVertex )
            {
                mVertexCoordinates.set_column( iVertex, mVertexCoordinates.get_column( tUsedVertexIndices( iVertex ) ) );

                // Update the vertex index map
                tVertexIndexMap[ tUsedVertexIndices( iVertex ) ] = iVertex;
            }

            // Trim the vertex coordinates matrix
            mVertexCoordinates.resize( mVertexCoordinates.n_rows(), tNumUsedVertices );

            // Loop through the facets and replace the vertex indices with the cleaned indices
            for ( auto& tFacet : mFacetConnectivity )
            {
                for ( moris_index& tVertexIndex : tFacet )
                {
                    tVertexIndex = tVertexIndexMap[ tVertexIndex ];
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::set_all_displacements( const Matrix< DDRMat >& aDisplacements )
    {
        MORIS_ASSERT( aDisplacements.n_rows() == this->get_spatial_dimension(), "Number of vertices in displacement matrix does not match number of vertices in mesh" );
        MORIS_ASSERT( aDisplacements.n_cols() == this->get_number_of_vertices(), "Number of dimensions in displacement matrix does not match number of dimensions in mesh" );

        mDisplacements = aDisplacements;

        // Update the normal vector for all the facets
        this->initialize_facet_normals();

#if MORIS_HAVE_ARBORX
        // Update the bounding volume hierarchy
        this->construct_bvh();
#endif
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

    const Matrix< DDRMat > Surface_Mesh::get_all_vertex_coordinates() const
    {
        return mVertexCoordinates + mDisplacements;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat > Surface_Mesh::get_vertex_coordinates( uint aVertexIndex ) const
    {
        return mVertexCoordinates.get_column( aVertexIndex ) + mDisplacements.get_column( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat > Surface_Mesh::get_all_original_vertex_coordinates() const
    {
        return mVertexCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat > Surface_Mesh::get_original_vertex_coordinates( const uint aVertexIndex ) const
    {
        return mVertexCoordinates.get_column( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Surface_Mesh::get_vertex_displacements() const
    {
        return mDisplacements;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Vector< moris_index > >& Surface_Mesh::get_facet_connectivity() const
    {
        return mFacetConnectivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Vector< moris_index > >& Surface_Mesh::get_vertex_connectivity() const
    {
        return mVertexConnectivity;
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

    real Surface_Mesh::get_intersection_tolerance() const
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
        // Assume a warning will be thrown, this will be changed if the raycast is successful
        bool tWarning = true;

        uint tDim = this->get_spatial_dimension();

        // Initialize random ray direction vector
        Matrix< DDRMat > tDirection( tDim, 1 );
        if ( tDim == 2 )
        {
            tDirection = { { 0.6398, -0.4472 } };
        }
        else
        {
            tDirection = { { 0.4990, -0.3534, -0.7912 } };
        }

        tDirection = tDirection / norm( tDirection );

        // Initialize vector to hold raycast result
        Intersection_Vector tIntersections;
        if ( tWarning )
        {
            // Initialize attempt counter
            uint tAttemptCounter = 0;

            // Cast the ray until it does not hit a warning or we reach a maximum number of attempts
            do
            {
                // Get a new random direction
                tDirection = this->random_direction();

                // Cast the ray
                tIntersections = this->cast_single_ray( aPoint, tDirection, tWarning, true );

                // Increment attempt counter
                tAttemptCounter++;

                // Check if we have exceeded the maximum number of attempts
                if ( tAttemptCounter > 100 )
                {
                    MORIS_LOG_WARNING( "Exceeded maximum number of attempts to resolve raycast warning. Setting to inferface" );
                    return Mesh_Region::INTERFACE;
                }
            } while ( tWarning );
        }

        // Determine the region based on the number of intersections or if any intersection is close to zero
        Mesh_Region tRegion = std::any_of( tIntersections.begin(), tIntersections.end(), [ this ]( std::pair< uint, real > aIntersection ) { return std::abs( aIntersection.second ) < mIntersectionTolerance; } )
                                    ? Mesh_Region::INTERFACE
                                    : static_cast< Mesh_Region >( tIntersections.size() % 2 );

        return tRegion;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Mesh_Region >
    Surface_Mesh::batch_get_region_from_raycast( Matrix< DDRMat >& aPoint ) const
    {
        // Get spatial dimension and number of points
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

        tDirection = tDirection / norm( tDirection );

        // Initialize output
        Vector< Mesh_Region > tRegions( tNumPoints, UNDEFINED );

        // Vector to track which rays have warnings
        Vector< Vector< bool > > tWarnings;

        // Preallocate for errors to minimize reallocations
        Matrix< DDRMat > tErroredOrigins( tDim, tNumPoints );    // Initially as large as possible
        Vector< uint >   tErroredIndices( tNumPoints );          // Tracks errored indices

        // Perform initial ray casting
        Vector< Vector< Intersection_Vector > > tIntersections = this->cast_batch_of_rays( aPoint, tDirection, tWarnings, false );

        // Filter the indices of errored points
        uint tNumWarnings = 0;
        for ( uint iRayIndex = 0; iRayIndex < tNumPoints; ++iRayIndex )
        {
            if ( tWarnings( iRayIndex )( 0 ) )
            {
                tErroredIndices( tNumWarnings ) = iRayIndex;
                tErroredOrigins.set_column( tNumWarnings++, aPoint.get_column( iRayIndex ) );
            }
        }

        // Recast all the errored rays until there are no more warnings
        if ( tNumWarnings > 0 )
        {
            // Initialize attempt counter
            uint tAttemptCounter = 0;

            do
            {
                MORIS_LOG_INFO( "%d rays failed to resolve", tNumWarnings );

                // Get a new random direction
                tDirection = this->random_direction();

                // Resize preallocated structures to match the actual number of errors
                tErroredOrigins.resize( tDim, tNumWarnings );

                // Create a temporary warnings vector for this batch of errored rays
                Vector< Vector< bool > > tNewWarnings;

                // Cast the errored rays
                Vector< Vector< Intersection_Vector > > tErroredIntersections = this->cast_batch_of_rays( tErroredOrigins, tDirection, tNewWarnings, false );

                // Update intersections for resolved rays and rebuild error list
                uint tNumNewWarnings = 0;
                for ( uint iWarning = 0; iWarning < tNumWarnings; ++iWarning )
                {
                    uint tRayIndex = tErroredIndices( iWarning );

                    if ( !tNewWarnings( iWarning )( 0 ) )    // No warning, means the ray was resolved
                    {
                        // Update intersection for resolved rays
                        tIntersections( tRayIndex )( 0 ) = tErroredIntersections( iWarning )( 0 );
                    }
                    else
                    {
                        // Keep unresolved rays in the error list
                        tErroredIndices( tNumNewWarnings ) = tRayIndex;
                        tErroredOrigins.set_column( tNumNewWarnings++, aPoint.get_column( tRayIndex ) );
                    }
                }

                // Update the number of warnings
                tNumWarnings = tNumNewWarnings;

                // Check if we have exceeded the maximum number of attempts
                if ( ++tAttemptCounter > 10 )
                {
                    MORIS_LOG_WARNING( "Exceeded maximum number of attempts to resolve raycast warnings. Setting these points to interface" );
                    for ( uint i = 0; i < tNumWarnings; ++i )
                    {
                        tRegions( tErroredIndices( i ) ) = Mesh_Region::INTERFACE;
                    }
                    break;
                }
            } while ( tNumWarnings > 0 );
        }

        // Loop through the intersections and determine the region for each point
        for ( uint iPoint = 0; iPoint < tNumPoints; iPoint++ )
        {
            // Check if this ray has already been classified
            if ( tRegions( iPoint ) == UNDEFINED )
            {
                // Get the intersections for this point FIXME: write a function that returns the distances only to avoid the std::pair overhead
                Intersection_Vector tIntersectionsForPoint = tIntersections( iPoint )( 0 );

                // Store the region for this point
                tRegions( iPoint ) = std::any_of( tIntersectionsForPoint.begin(), tIntersectionsForPoint.end(), [ this ]( std::pair< uint, real > aIntersection ) { return std::abs( aIntersection.second ) < mIntersectionTolerance; } )
                                           ? Mesh_Region::INTERFACE
                                           : static_cast< Mesh_Region >( tIntersectionsForPoint.size() % 2 );
            }
        }

        return tRegions;
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Vector
    Surface_Mesh::cast_single_ray(
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection,
            bool&                   aWarning,
            bool                    aIgnoreWarnings ) const
    {
#if MORIS_HAVE_ARBORX
        // Get the facets that the ray could intersect using arborx
        Vector< uint > tCandidateFacets = this->preselect_with_arborx( aPoint, aDirection );
#else
        // Get all of the facets
        Vector< uint > tCandidateFacets( this->get_number_of_facets() );
        std::iota( tCandidateFacets.begin(), tCandidateFacets.end(), 0 );
#endif

        // Compute the intersection locations for all of the candidates, remove duplicates, and sort
        return this->determine_valid_intersections_from_candidates( aPoint, aDirection, tCandidateFacets, aWarning, aIgnoreWarnings );
    }


    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< Intersection_Vector > >
    Surface_Mesh::cast_batch_of_rays(
            Matrix< DDRMat >&         aOrigins,
            Matrix< DDRMat >&         aDirections,
            Vector< Vector< bool > >& aWarnings,
            bool                      aIgnoreWarnings ) const
    {
        // Get the number of origins
        const uint tNumberOfOrigins    = aOrigins.n_cols();
        const uint tNumberOfDirections = aDirections.n_cols();

        // Prepare warning vector if needed
        if ( not aIgnoreWarnings )
        {
            aWarnings.resize( tNumberOfOrigins, Vector< bool >( tNumberOfDirections, false ) );
        }

#if MORIS_HAVE_ARBORX
        // Get the facets that the ray could intersect
        Vector< Vector< Vector< uint > > > tCandidateFacets = this->batch_preselect_with_arborx( aOrigins, aDirections );
#else

        // initialize the candidate vector
        Vector< Vector< Vector< uint > > > tCandidateFacets( tNumberOfOrigins, Vector< Vector< uint > >( tNumberOfDirections ) );

        // Get all of the facets
        Vector< uint > tAllFacetIndices( this->get_number_of_facets() );
        std::iota( tAllFacetIndices.begin(), tAllFacetIndices.end(), 0 );

        // Assign all of the facets to the candidate vector
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            for ( uint iDirection = 0; iDirection < tNumberOfDirections; iDirection++ )
            {
                tCandidateFacets( iOrigin )( iDirection ) = tAllFacetIndices;
            }
        }
#endif

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
                bool tWarning;

                // Get the candidate facets for this ray
                Vector< uint > tCandidateFacetsForRay = tCandidateFacets( iOrigin )( iDirection );

                // Compute which facets are intersected, remove duplicates, sort, and store in output
                tAllIntersections( iOrigin )( iDirection ) = this->determine_valid_intersections_from_candidates(
                        aOrigins.get_column( iOrigin ),
                        aDirections.get_column( iDirection ),
                        tCandidateFacetsForRay,
                        tWarning,
                        aIgnoreWarnings );

                // Store the warning if needed
                if ( not aIgnoreWarnings )
                {
                    aWarnings( iOrigin )( iDirection ) = tWarning;
                }
            }
        }

        return tAllIntersections;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< Intersection_Vector > >
    Surface_Mesh::cast_batch_of_rays(
            Matrix< DDRMat >&           aOrigins,
            Vector< Matrix< DDRMat > >& aDirections,
            Vector< Vector< bool > >&   aWarnings,
            bool                        aIgnoreWarnings ) const
    {
        // Get the number of origins
        uint tNumberOfOrigins = aOrigins.n_cols();

        MORIS_ASSERT( tNumberOfOrigins == aDirections.size(), "To cast a batch of rays with different directions for each ray, the size of the vector of directions (%lu) must match the columns of aOrigins (%d)", aDirections.size(), tNumberOfOrigins );

        // Prepare warning vector if needed
        if ( not aIgnoreWarnings )
        {
            aWarnings.resize( tNumberOfOrigins, false );
        }

#if MORIS_HAVE_ARBORX
        // Get the facets that the ray could intersect
        Vector< Vector< Vector< uint > > > tCandidateFacets = this->batch_preselect_with_arborx( aOrigins, aDirections );

        // Resize the warnings vector
        if ( not aIgnoreWarnings )
        {
            for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
            {
                aWarnings( iOrigin ).resize( aDirections( iOrigin ).n_cols() );
            }
        }
#else

        // initialize the candidate vector
        Vector< Vector< Vector< uint > > > tCandidateFacets( tNumberOfOrigins );

        // Get all of the facets
        Vector< uint > tAllFacetIndices( this->get_number_of_facets() );
        std::iota( tAllFacetIndices.begin(), tAllFacetIndices.end(), 0 );

        // Assign all of the facets to the candidate vector
        for ( uint iOrigin = 0; iOrigin < tNumberOfOrigins; iOrigin++ )
        {
            // Get the number of directions for this origin
            uint tNumberOfDirections = aDirections( iOrigin ).n_cols();

            // resize the candidates and the warnings
            tCandidateFacets( iOrigin ).resize( tNumberOfDirections );
            if ( not aIgnoreWarnings )
            {
                aWarnings( iOrigin ).resize( tNumberOfDirections );
            }

            // fill the candidate vector with all facets
            for ( uint iDirection = 0; iDirection < tNumberOfDirections; iDirection++ )
            {
                tCandidateFacets( iOrigin )( iDirection ) = tAllFacetIndices;
            }
        }
#endif

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
                bool tWarning;

                // Get the candidate facets for this ray
                Vector< uint > tCandidateFacetsForRay = tCandidateFacets( iOrigin )( iDirection );

                // Compute which facets are intersected, remove duplicates, sort, and store in output
                tAllIntersections( iOrigin )( iDirection ) = this->determine_valid_intersections_from_candidates(
                        aOrigins.get_column( iOrigin ),
                        aDirections( iOrigin ).get_column( iDirection ),
                        tCandidateFacetsForRay,
                        tWarning,
                        aIgnoreWarnings );

                // Store the warning if needed
                if ( not aIgnoreWarnings )
                {
                    aWarnings( iOrigin )( iDirection ) = tWarning;
                }
            }
        }

        return tAllIntersections;
    }

    // --------------------------------------------------------------------------------------------------------------

    Intersection_Vector Surface_Mesh::determine_valid_intersections_from_candidates(
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection,
            const Vector< uint >&   aCandidateFacets,
            bool&                   aWarning,
            bool                    aIgnoreWarnings ) const
    {
        // Assume this ray will not be a pathological case
        aWarning = false;

        Intersection_Vector tIntersections( aCandidateFacets.size() );
        uint                tNumberOfValidIntersections = 0;

        for ( uint iCandidate : aCandidateFacets )
        {
            // Compute the intersection location
            real tIntersection = this->moller_trumbore( iCandidate, aPoint, aDirection, aWarning );

            // If we are not ignoring warnings, check for one
            if ( not aIgnoreWarnings and aWarning )
            {
                // break the rest of the postprocessing
                tNumberOfValidIntersections = 0;
                aWarning                    = true;
                break;
            }

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

            return this->sort_and_find_unique_intersections( tIntersections );
        }
        else
        {
            return {};
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    real Surface_Mesh::moller_trumbore(
            uint                    aFacet,
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection,
            bool&                   aWarning ) const
    {
        switch ( this->get_spatial_dimension() )
        {
            case 2:
            {
                return moller_trumbore_2D( aFacet, aPoint, aDirection, aWarning );
            }
            case 3:
            {
                return moller_trumbore_3D( aFacet, aPoint, aDirection, aWarning );
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
            const Matrix< DDRMat >& aDirection,
            bool&                   aWarning ) const
    {
        // Assume this ray will not hit a vertex
        aWarning = false;

        // Get the vertex coordinates for the requested facet
        Matrix< DDRMat > tFacetCoordinates = this->get_all_vertex_coordinates_of_facet( aFacet );

        // Build the vector pointing along the facet edge
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

        // Get the determinant of the edge and the cast direction
        real tDet = cross_2d( aDirection, tEdge );

        // If the determinant is close to zero, the ray is parallel or colinear to the line
        if ( std::abs( tDet ) < mIntersectionTolerance )
        {
            // Check for colinearity
            if ( std::abs( cross_2d( tRayToVertex, aDirection ) < mIntersectionTolerance ) )
            {
                // Check if the point is in between the vertices
                if ( ( tFacetCoordinates( 0, 0 ) - aPoint( 0 ) ) * ( tFacetCoordinates( 0, 1 ) - aPoint( 0 ) ) < 0.0
                        or ( tFacetCoordinates( 1, 0 ) - aPoint( 1 ) ) * ( tFacetCoordinates( 1, 1 ) - aPoint( 1 ) ) < 0.0 )
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
        real tInverseDet = 1.0 / tDet;

        // Solve the 2D system
        real tDistance = cross_2d( tRayToVertex, tEdge ) * tInverseDet;
        real tU        = cross_2d( tRayToVertex, aDirection ) * tInverseDet;

        // Check the u parameter for edge cases
        aWarning = ( ( std::abs( tU ) < 1e-3 and std::abs( tU ) > mIntersectionTolerance )
                           or ( std::abs( tU - 1.0 ) < 1e-3 and std::abs( tU - 1.0 ) > mIntersectionTolerance ) )
               and std::abs( tDistance ) > mIntersectionTolerance;


        // Check if the intersection is within the line segment
        if ( tU < -mIntersectionTolerance or tU > 1.0 + mIntersectionTolerance or tDistance < -mIntersectionTolerance )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }

        return ( std::abs( tDistance ) < mIntersectionTolerance ? 0.0 : tDistance );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh::moller_trumbore_3D(
            uint                    aFacet,
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aDirection,
            bool&                   aWarning ) const
    {
        // Assume this ray will not hit a vertex
        aWarning = false;

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
        Matrix< DDRMat > tRayToVertex;
        if ( aPoint.n_cols() == 1 )
        {
            tRayToVertex = aPoint - tVertexCoordinates.get_column( 0 );
        }
        else
        {
            tRayToVertex = trans( aPoint ) - tVertexCoordinates.get_column( 0 );
        }

        // Compute the u parameter
        real tU = dot( tRayToVertex, tP ) * tInverseDeterminant;

        // If the u parameter is < 0.0 or > 1.0, the intersection is outside the triangle
        if ( tU < -mIntersectionTolerance or tU > 1.0 + mIntersectionTolerance )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }

        // Compute the vector from the origin to the second vertex
        Matrix< DDRMat > tQ = cross( tRayToVertex, tEdge1 );

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
        if ( tDistance < -mIntersectionTolerance )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }
        else if ( tDistance < mIntersectionTolerance )
        {
            return 0.0;
        }
        else
        {
            // Check if the u or v parameters are close to an edge case and throw a warning if so
            aWarning = ( ( std::abs( tU ) < 1e-3 and std::abs( tU ) > mIntersectionTolerance ) or ( std::abs( tU - 1.0 ) < 1e-3 and std::abs( tU - 1.0 ) > mIntersectionTolerance ) )
                    or ( ( std::abs( tV ) < 1e-3 and std::abs( tV ) > mIntersectionTolerance ) or ( std::abs( tV - 1.0 ) < 1e-3 and std::abs( tV - 1.0 ) > mIntersectionTolerance ) );

            return tDistance;
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Intersection_Vector Surface_Mesh::sort_and_find_unique_intersections( Intersection_Vector& aIntersections ) const
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

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::random_direction() const
    {
        // Random device for random number seeding
        std::random_device tRD;

        // Pseudo-random number generator (Mersenne Twister)
        std::mt19937 tGen( tRD() );

        // Define a uniform real distribution in the range [-1.0, 1.0)
        std::uniform_real_distribution< double > tDist( -1.0, 1.0 );

        uint tDim = this->get_spatial_dimension();

        // Get a random number for each direction
        Matrix< DDRMat > tDirection( tDim, 1 );
        for ( uint iDim = 0; iDim < tDim; iDim++ )
        {
            tDirection( iDim ) = tDist( tGen );
        }

        // Normalize the direction
        return tDirection / norm( tDirection );
    }

    //--------------------------------------------------------------------------------------------------------------


#if MORIS_HAVE_ARBORX
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
                        const auto tOrigin = arborx::coordinate_to_arborx_point< ArborX::Point >( aOrigins.get_column( iOrigin ) );

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
                    const auto tOrigin = arborx::coordinate_to_arborx_point< ArborX::Point >( aOrigins.get_column( iOrigin ) );

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

                // Loop through the directions for this origin point
                for ( uint iDirection = 0; iDirection < tNumDirections; iDirection++ )
                {
                    Vector< uint > tFacetIndicesForRay( tOffsets( tOffsetIndex + 1 ) - tOffsets( tOffsetIndex ) );

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
                // Loop through the directions for this origin point
                for ( uint iDirection = 0; iDirection < tNumDirections; iDirection++ )
                {
                    Vector< uint > tFacetIndicesForRay( tOffsets( tOffsetIndex + 1 ) - tOffsets( tOffsetIndex ) );

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

#endif

    // --------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::write_to_file( const std::string& aFilePath ) const
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
            tFile << "\n";
        }

        for ( uint iFacet = 0; iFacet < this->get_number_of_facets(); iFacet++ )
        {
            tFile << "f ";
            const Vector< moris_index > tVertexIndices = this->get_facets_vertex_indices( iFacet );
            for ( uint iVertexIndex = 0; iVertexIndex < tVertexIndices.size(); iVertexIndex++ )
            {
                tFile << tVertexIndices( iVertexIndex ) + 1 << " ";
            }
            tFile << "\n";
        }

        // close file
        tFile.close();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::build_vertex_connectivity()
    {
        // Initialize a set to store unique vertex connectivity for every vertex
        uint                              tNumVertices = this->get_number_of_vertices();
        Vector< std::set< moris_index > > tVertexConnectivity( tNumVertices );

        uint tNumVertsPerFacet = this->get_facets_vertex_indices( 0 ).size();

        // Loop through all facets of the surface mesh
        for ( auto& tFacet : this->get_facet_connectivity() )
        {
            // Loop through all vertices of the facet
            for ( uint iVertex = 0; iVertex < tNumVertsPerFacet; iVertex++ )
            {
                // Get the current vertex and the next vertex (wrapping around)
                moris_index tCurrentVertex = tFacet( iVertex );
                moris_index tNextVertex    = tFacet( ( iVertex + 1 ) % tNumVertsPerFacet );

                // Add both to each others sets
                tVertexConnectivity( tCurrentVertex ).insert( tNextVertex );
                tVertexConnectivity( tNextVertex ).insert( tCurrentVertex );
            }
        }

        // Convert the sets to vectors and store them in the vertex connectivity
        mVertexConnectivity.resize( tNumVertices );
        for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
        {
            // Convert the set to a vector
            mVertexConnectivity( iVertex ) = Vector< moris_index >( tVertexConnectivity( iVertex ).begin(), tVertexConnectivity( iVertex ).end() );
        }
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

    // --------------------------------------------------------------------------------------------------------------
    // Quantities of Interest
    // --------------------------------------------------------------------------------------------------------------

    real Surface_Mesh::compute_volume() const
    {
        real tVolume = 0.0;

        // Apply shoelace formula for 2D surface mesh BRENDAN NOT GENERALIZABLE
        if ( this->get_spatial_dimension() == 2 )
        {
            for ( uint iFacet = 0; iFacet < this->get_number_of_facets(); iFacet++ )
            {
                Matrix< DDRMat > tCoords = this->get_all_vertex_coordinates_of_facet( iFacet );

                tVolume += 0.5 * ( tCoords( 0, 0 ) * tCoords( 1, 1 ) - tCoords( 0, 1 ) * tCoords( 1, 0 ) );
            }
        }
        else
        {
            MORIS_ERROR( false, "Surface Mesh volume computation only implemented for 2D (lines) meshes" );
        }

        return tVolume;
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_dvolume_dvertex( uint aVertexIndex ) const
    {
        uint tDim = this->get_spatial_dimension();

        Matrix< DDRMat > tDVolumeDVertex( 1, tDim, MORIS_REAL_MAX );

        if ( this->get_spatial_dimension() == 2 )
        {
            // Get the vertices connected to this vertex
            const Vector< moris_index >& tNeighbors = mVertexConnectivity( aVertexIndex );

            // Determine the order of the vertices
            moris_index tLowNeighbor  = tNeighbors( 0 ) < tNeighbors( 1 ) ? tNeighbors( 0 ) : tNeighbors( 1 );
            moris_index tHighNeighbor = tNeighbors( 0 ) > tNeighbors( 1 ) ? tNeighbors( 0 ) : tNeighbors( 1 );

            // Get vertex coordinates
            const Matrix< DDRMat > tLowCoords  = this->get_vertex_coordinates( tLowNeighbor );
            const Matrix< DDRMat > tHighCoords = this->get_vertex_coordinates( tHighNeighbor );

            tDVolumeDVertex( 0 ) = 0.5 * ( tHighCoords( 1 ) - tLowCoords( 1 ) );
            tDVolumeDVertex( 1 ) = 0.5 * ( tLowCoords( 0 ) - tHighCoords( 0 ) );
        }
        else
        {
            MORIS_ERROR( false, "Surface Mesh dVolume/dVertex computation only implemented for 2D (lines) meshes" );
        }

        return tDVolumeDVertex;
    }


}    // namespace moris::mtk
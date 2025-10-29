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

#include "SDF_Tools.hpp"

#include <set>
#include <random>
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"
#include "fn_trans.hpp"
#include "fn_eye.hpp"
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
            , mFacetToVertexConnectivity( aFacetConnnectivity )
            , mIntersectionTolerance( aIntersectionTolerance )
    {
        // Remove the unused vertices, update the facet connectivity
        this->clean_extraneous_vertices();

        // Initialize distortion vectors/matrices
        this->reset_coordinates();

        // Compute the normals of the facets
        this->initialize_facet_normals();

        // Build the vertex to vertex connectivity
        this->build_vertex_connectivity();

        // Build the vertex to facet connectivity
        this->build_vertex_to_facet_connectivity();

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
        // NOTE: relies on surface mesh being lines in 2D and triangles in 3D
        Vector< moris_index > tUsedVertexIndices( mFacetToVertexConnectivity.size() * this->get_spatial_dimension(), MORIS_INDEX_MAX );

        // Loop through all the facets and store the vertices they use
        uint tIndex = 0;
        for ( auto tFacet : mFacetToVertexConnectivity )
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
            for ( auto& tFacet : mFacetToVertexConnectivity )
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

        // Recompute the normals of the facets connected to this vertex
        for ( moris_index iF : this->get_vertexs_facet_indices( aVertexIndex ) )
        {
            this->compute_facet_normal( iF );
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
        return mFacetToVertexConnectivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Vector< moris_index > >& Surface_Mesh::get_vertex_connectivity() const
    {
        return mVertexToVertexConnectivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< moris_index >& Surface_Mesh::get_facets_vertex_indices( const uint aFacetIndex ) const
    {
        MORIS_ASSERT( aFacetIndex < this->get_number_of_facets(), "Surface_Mesh::get_facets_vertex_indices() - Facet index %d out of bounds (mesh has %d facets)", aFacetIndex, this->get_number_of_facets() );

        return mFacetToVertexConnectivity( aFacetIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< moris_index >& Surface_Mesh::get_vertexs_facet_indices( const uint aVertexIndex ) const
    {
        MORIS_ASSERT( aVertexIndex < this->get_number_of_vertices(), "Surface_Mesh::get_vertexs_facet_indices() - Vertex index %d out of bounds (mesh has %d vertices)", aVertexIndex, this->get_number_of_vertices() );

        return mVertexToFacetConnectivity( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< moris_index >& Surface_Mesh::get_vertex_neighbors( const uint aVertexIndex ) const
    {
        MORIS_ASSERT( aVertexIndex < this->get_number_of_vertices(), "Surface_Mesh::get_vertex_neighbors() - Vertex index %d out of bounds (mesh has %d vertices)", aVertexIndex, this->get_number_of_vertices() );

        return mVertexToVertexConnectivity( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::get_all_vertex_coordinates_of_facet( uint aFacetIndex ) const
    {
        // Get the facets vertex indices
        Vector< moris_index > tVertices    = this->get_facets_vertex_indices( aFacetIndex );
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
        return mFacetToVertexConnectivity.size();
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
        Matrix< DDRMat > tOriginToV0;
        if ( aPoint.n_cols() == 1 )
        {
            tOriginToV0 = tFacetCoordinates.get_column( 0 ) - aPoint;
        }
        else
        {
            tOriginToV0 = tFacetCoordinates.get_column( 0 ) - trans( aPoint );
        }

        // Get the determinant of the edge and the cast direction
        real tDet = cross_2d( aDirection, tEdge );

        // If the determinant is close to zero, the ray is parallel or colinear to the line
        if ( std::abs( tDet ) < mIntersectionTolerance )
        {
            // Check for colinearity
            if ( std::abs( cross_2d( tOriginToV0, aDirection ) < mIntersectionTolerance ) )
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
        real tDistance = cross_2d( tOriginToV0, tEdge ) * tInverseDet;
        real tU        = cross_2d( tOriginToV0, aDirection ) * tInverseDet;

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
        Matrix< DDRMat > tOriginToV0;
        if ( aPoint.n_cols() == 1 )
        {
            tOriginToV0 = aPoint - tVertexCoordinates.get_column( 0 );
        }
        else
        {
            tOriginToV0 = trans( aPoint ) - tVertexCoordinates.get_column( 0 );
        }

        // Compute the u parameter
        real tU = dot( tOriginToV0, tP ) * tInverseDeterminant;

        // If the u parameter is < 0.0 or > 1.0, the intersection is outside the triangle
        if ( tU < -mIntersectionTolerance or tU > 1.0 + mIntersectionTolerance )
        {
            return std::numeric_limits< real >::quiet_NaN();
        }

        // Compute the vector from the origin to the second vertex
        Matrix< DDRMat > tQ = cross( tOriginToV0, tEdge1 );

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
        mVertexToVertexConnectivity.resize( tNumVertices );
        for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
        {
            // Convert the set to a vector
            mVertexToVertexConnectivity( iVertex ) = Vector< moris_index >( tVertexConnectivity( iVertex ).begin(), tVertexConnectivity( iVertex ).end() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::build_vertex_to_facet_connectivity()
    {
        mVertexToFacetConnectivity.resize( this->get_number_of_vertices(), Vector< moris_index >( this->get_spatial_dimension() ) );

        for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
        {
            Vector< moris_index > tFacetVertices = this->get_facets_vertex_indices( iF );

            for ( uint iV = 0; iV < tFacetVertices.size(); iV++ )
            {
                mVertexToFacetConnectivity( tFacetVertices( iV ) )( iV ) = (moris_index)iF;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::set_vertex_coordinates( const uint aVertexIndex, const Matrix< DDRMat >& aCoordinates )
    {
        mVertexCoordinates.set_column( aVertexIndex, aCoordinates );
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_facet_normal( const uint aFacetIndex )
    {
        uint const tDim = this->get_spatial_dimension();

        Matrix< DDRMat > tNormal( tDim, 1 );
        switch ( tDim )
        {
            case 2:
            {
                Matrix< DDRMat > tCoords = this->get_all_vertex_coordinates_of_facet( aFacetIndex );

                // { { tY2 - tY1  }, { tX1 - tX2 } }
                tNormal( 0 ) = tCoords( 1, 1 ) - tCoords( 1, 0 );
                tNormal( 1 ) = tCoords( 0, 0 ) - tCoords( 0, 1 );
                tNormal      = tNormal / norm( tNormal );
                break;
            }
            case 3:
            {
                Matrix< DDRMat > tCoords = this->get_all_vertex_coordinates_of_facet( aFacetIndex );

                tNormal = cross( tCoords.get_column( 1 ) - tCoords.get_column( 0 ), tCoords.get_column( 2 ) - tCoords.get_column( 0 ) );
                tNormal = tNormal / norm( tNormal );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Surface Mesh facet normals only implemented for 2D (lines) or 3D (triangles) meshes" );
            }
        }

        mFacetNormals.set_column( aFacetIndex, tNormal );

        return tNormal;
    }

    // -------------------------------------------------------------------------------------------------------------

    void Surface_Mesh::initialize_facet_normals()
    {
        auto const tNumFacets = static_cast< moris::size_t >( mFacetToVertexConnectivity.size() );
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

    Matrix< DDRMat > Surface_Mesh::compute_facet_centroid( const uint aFacetIndex ) const
    {
        uint             tDim = this->get_spatial_dimension();
        Matrix< DDRMat > tCentroid( tDim, 1 );

        // Get vertex coordinates of the facet
        Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( aFacetIndex );

        // Loop over vertices in the facet
        for ( uint iV = 0; iV < tVertexCoordinates.n_cols(); iV++ )
        {
            // Loop over dimensions to accumulate coordinates
            for ( uint iDim = 0; iDim < tDim; iDim++ )
            {
                // Compute centroid as the average of vertex coordinates
                tCentroid( iDim ) += tVertexCoordinates( iDim, iV );
            }
        }

        // Divide by number of vertices to get the average
        // WARNING: This assumes all facets are lines in 2D or triangles in 3D
        return tCentroid / (real)tVertexCoordinates.n_cols();
    }

    // -------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_facet_centroids() const
    {
        // Get number of facets and spatial dimension
        uint tNumFacets = this->get_number_of_facets();
        uint tDim       = this->get_spatial_dimension();

        // Initialize matrix to hold centroids
        Matrix< DDRMat > tCentroids( tDim, tNumFacets );

        // Loop over facets to compute centroids
        for ( uint iF = 0; iF < tNumFacets; iF++ )
        {
            tCentroids.set_column( iF, this->compute_facet_centroid( iF ) );
        }

        return tCentroids;
    }

    // --------------------------------------------------------------------------------------------------------------

    Ray_Cones Surface_Mesh::build_ray_cone_angles( real aConeAngle, uint aNumPolarRays, uint aNumAzimuthRays ) const
    {
        const uint tNumFacets = this->get_number_of_facets();
        const uint tDim       = this->get_spatial_dimension();

        // Initialize return variable
        Ray_Cones tRayCones( tNumFacets, aNumPolarRays * aNumAzimuthRays, tDim );

        switch ( tDim )
        {
            case 2:
            {
                const real tdAlpha = aConeAngle / static_cast< real >( aNumPolarRays );

                Matrix< DDRMat > tRotation;

                for ( uint iF = 0; iF < tNumFacets; iF++ )
                {
                    Matrix< DDRMat > tFacetNormal = this->get_facet_normal( iF );

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

                        tRayCones.mTheta( iR )            = tRayAngle * M_PI / 180;    // Convert to radians
                        tRayCones.mDirectionWeights( iR ) = 1.0 / std::abs( tRayAngle );

                        tRotation = sdf::rotation_matrix( tRayCones.mTheta( iR ) );

                        // The surface mesh has outward normals. Thus, we need to invert the cone direction to shoot inward
                        Matrix< DDRMat > tConeRay = -tRotation * tFacetNormal;

                        // Check that the ray is opposite of the facet normal
                        MORIS_ASSERT( dot( tConeRay, tFacetNormal ) < 0, "Surface_Mesh::build_ray_cone_angles - Ray direction for cone is not opposite to facet normal." );

                        // Check that the norm is not too small
                        MORIS_ERROR( norm( tConeRay ) >= 1e-9, "Surface_Mesh::build_ray_cone_angles - Ray direction has norm <1e-9. Should be unit." );

                        // if ( dot( tConeRay, tFacetNormal ) > 0 )
                        // {
                        //     tConeRay = -1.0 * tConeRay;
                        // }

                        // if ( norm( tConeRay ) < 1e-9 )
                        // {
                        //     for ( uint iMatSize = 0; iMatSize < 2; iMatSize++ )
                        //     {
                        //         tConeRay( iMatSize, 0 ) += 1e-12;
                        //     }
                        // }

                        // // Apply the tolerance to avoid -0.0
                        // for ( uint iConeRay = 0; iConeRay < tConeRay.n_rows(); iConeRay++ )
                        // {
                        //     if ( std::abs( tConeRay( iConeRay, 0 ) ) < 1e-9 )
                        //     {
                        //         tConeRay( iConeRay, 0 ) = 0.000000000000000e+00;    // Set to exact zero if close to zero
                        //     }
                        // }

                        MORIS_ASSERT( norm( tConeRay ) - 1.0 < 1e-12, "Surface_Mesh::build_ray_cone_angles - Ray direction for cone is not unit." );

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
                    Matrix< DDRMat > e3 = this->get_facet_normal( iF );    // Facet normal is the e3 vector

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

        // // brendan delete
        // std::cout << "Ray normals:\n";
        // std::cout << "[ ";
        // for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
        // {
        //     for ( uint iD = 0; iD < 2; iD++ )
        //     {
        //         std::cout << this->get_facet_normal( iF )( iD ) << " ";
        //     }
        //     std::cout << "; ";
        // }
        // std::cout << "]\n";

        // std::cout << "Ray directions:\n";
        // std::cout << "{ ";
        // for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
        // {
        //     std::cout << "[ ";
        //     for ( uint iR = 0; iR < aNumPolarRays * aNumAzimuthRays; iR++ )
        //     {
        //         for ( uint iD = 0; iD < 2; iD++ )
        //         {
        //             std::cout << tRayCones.mRayDirections( iF )( iD, iR ) << " ";
        //         }
        //         std::cout << "; ";
        //     }
        //     std::cout << "], ";
        // }
        // std::cout << "}\n";

        // this->write_to_file( "integration_surface_mesh.obj" );

        return tRayCones;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_cone_direction_from_angle( const Matrix< DDRMat >& aNormal, real aTheta, real aPhi ) const
    {
        uint tDim = this->get_spatial_dimension();

        Matrix< DDRMat > tDirection( tDim, 1 );

        switch ( tDim )
        {
            case 2:
            {
                Matrix< DDRMat > tRotation = sdf::rotation_matrix( aTheta );    // Already in radians

                tDirection = -tRotation * aNormal;

                break;
            }
            case 3:
            {
                MORIS_ERROR( false, "Not implemented yet for 3D." );    // brendan FIXME

                break;
            }
            default:
                MORIS_ERROR( false, "Only 2D-3D implementation" );
                break;
        }

        return tDirection;
    }

    // --------------------------------------------------------------------------------------------------------------

    Vector< Intersection_Vector > Surface_Mesh::determine_nearest_nontrivial_intersections( Vector< Vector< Intersection_Vector > >& aAllIntersections ) const
    {
        MORIS_ASSERT( aAllIntersections.size() == this->get_number_of_facets(), "Input intersections size does not match number of facets in the mesh." );

        // Initialize return variable
        Vector< Intersection_Vector > tNearestIntersections( this->get_number_of_facets(), Intersection_Vector( aAllIntersections( 0 ).size() ) );

        // Loop over facets
        for ( uint iF = 0; iF < aAllIntersections.size(); iF++ )
        {
            MORIS_ASSERT( aAllIntersections( iF ).size() == aAllIntersections( 0 ).size(), "Inconsistent number of rays per facet in the input intersections." );

            // Loop over rays for this facet
            for ( uint iR = 0; iR < aAllIntersections( iF ).size(); iR++ )
            {
                MORIS_ASSERT( aAllIntersections( iF )( iR ).size() > 0, "No intersections found for facet %d, ray %d. Is the normal correct?", iF, iR );

                // Loop over the intersections and find the first intersection that is not with the facet the ray came from
                for ( uint iI = 0; iI < aAllIntersections( iF )( iR ).size(); iI++ )
                {
                    if ( aAllIntersections( iF )( iR )( iI ).first != iF )    // Check that we didn't get an intersection with the originating facet by mistake
                    {
                        tNearestIntersections( iF )( iR ) = aAllIntersections( iF )( iR )( iI );
                        break;
                    }

                    MORIS_ASSERT( iI < aAllIntersections( iF )( iR ).size() - 1, "No non-trivial intersection found for facet %d, ray %d.", iF, iR );
                }
            }
        }

        return tNearestIntersections;
    }

    // --------------------------------------------------------------------------------------------------------------

    Shape_Diameter_Distances Surface_Mesh::cast_shape_diameter_ray_cones( real aConeAngle, uint aNumPolarRays, uint aNumAzimuthRays ) const
    {
        // Build a cone of rays for each vertex centered around its vertex normal
        Ray_Cones tRayCones = this->build_ray_cone_angles( aConeAngle, aNumPolarRays, aNumAzimuthRays );

        // Get the current vertex coordinates for the entire surface mesh
        Matrix< DDRMat > tRayOrigins = this->compute_facet_centroids();

        // Initialize the output struct - get the direction weights from the ray cones
        Shape_Diameter_Distances tConeDistances( tRayCones, this->get_number_of_vertices(), aNumPolarRays * aNumAzimuthRays );

        // Batch process all rays for all vertices
        Vector< Vector< bool > >                tWarnings;
        Vector< Vector< Intersection_Vector > > tIntersections = this->cast_batch_of_rays( tRayOrigins, tRayCones.mRayDirections, tWarnings );

        // Get just the nearest intersection for each ray
        tConeDistances.mDistances = this->determine_nearest_nontrivial_intersections( tIntersections );

        // // brendan delete
        // std::cout << "Ray distances:\n{ ";
        // for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
        // {
        //     std::cout << "[ ";
        //     for ( uint iR = 0; iR < aNumPolarRays * aNumAzimuthRays; iR++ )
        //     {
        //         std::cout << tConeDistances.mDistances( iF )( iR ).second << " ";
        //     }
        //     std::cout << "], ";
        // }
        // std::cout << "};\n";

        // std::cout << "Ray facets:\n{ ";
        // for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
        // {
        //     std::cout << "[ ";
        //     for ( uint iR = 0; iR < aNumPolarRays * aNumAzimuthRays; iR++ )
        //     {
        //         std::cout << tConeDistances.mDistances( iF )( iR ).first << " ";
        //     }
        //     std::cout << "], ";
        // }
        // std::cout << "};\n";

        // std::cout << "weights:\n[ ";
        // for ( uint iR = 0; iR < aNumPolarRays * aNumAzimuthRays; iR++ )
        // {
        //     std::cout << tConeDistances.mRayCones.mDirectionWeights( iR ) << ", ";
        // }
        // std::cout << "];\n";

        return tConeDistances;
    }

    //-------------------------------------------------------------------------------------------------------------

    Vector< real > Surface_Mesh::compute_nodal_shape_diameter( real aConeAngle, uint aNumPolarRays, uint aNumAzimuthRays )
    {
        // Get the minimum distance for every ray cone at every facet
        Shape_Diameter_Distances tConeDistances = this->cast_shape_diameter_ray_cones( aConeAngle, aNumPolarRays, aNumAzimuthRays );

        // Get the number of facets in the mesh
        uint tNumFacets = this->get_number_of_facets();

        // Get the number of rays per facet
        uint tNumRays = tConeDistances.mRayCones.mDirectionWeights.size();

        // Check that we got a distance for every facet and ray
        MORIS_ASSERT( tConeDistances.mDistances.size() == tNumFacets, "Inconsistent number of facets in the shape diameter distances." );

        // Centroids needed for raycast sensitivities
        Matrix< DDRMat > tCentroids = this->compute_facet_centroids();

        // Initialize return variable
        Vector< real > tShapeDiameterValues( tNumFacets, MORIS_REAL_MAX );

        // Reset the derivative to store for later
        mdShapeDiameterdVertex.set_size( this->get_number_of_vertices(), this->get_spatial_dimension(), 0.0 );

        // Get the facet measures for sensitivity calculations
        Vector< real > tFacetMeasures = this->compute_facet_measure();

        // Loop through every facet in this surface mesh
        // std::cout << "[ ";    // brendan delete prints
        for ( uint iF = 0; iF < tNumFacets; iF++ )
        {
            MORIS_ASSERT( tConeDistances.mDistances( iF ).size() == tNumRays, "Inconsistent number of rays in cone for facet %d (Expected %d, Got %ld)", iF, tNumRays, tConeDistances.mDistances( iF ).size() );

            // Get the ray distances for this facet (one distance per ray)
            Intersection_Vector& tRaysOnFacet = tConeDistances.mDistances( iF );

            Intersection_Vector tSortedRaysOnFacet = tRaysOnFacet;    // FIXME annoying copy

            // Compute the median ray distance for this facet
            std::sort( tSortedRaysOnFacet.begin(), tSortedRaysOnFacet.end(), []( const std::pair< moris_index, real >& aR1, const std::pair< moris_index, real >& aR2 ) { return aR1.second < aR2.second; } );
            real tMedian;
            if ( tNumRays % 2 == 0 )
            {
                tMedian = ( tSortedRaysOnFacet( tNumRays / 2 - 1 ).second + tSortedRaysOnFacet( tNumRays / 2 ).second ) / 2;
            }
            else
            {
                tMedian = tSortedRaysOnFacet( tNumRays / 2 ).second;
            }

            // Compute the mean ray distance for this facet
            real tMean = 0.0;
            for ( const auto& tDistance : tRaysOnFacet )
            {
                tMean += tDistance.second;
            }
            tMean = tMean / tNumRays;

            // Compute the standard deviation ray distances for this facet
            real tStdDeviation = 0.0;
            for ( const auto& tDistance : tRaysOnFacet )
            {
                tStdDeviation += ( tDistance.second - tMean ) * ( tDistance.second - tMean );    // variance
            }
            tStdDeviation = std::sqrt( tStdDeviation / (real)tNumRays );    // take sqrt for std deviation

            // Compute the weighted sum for entries within one standard deviation from the median
            real tDiamWeightedSum   = 0.0;
            real tSumOfValidWeights = 0.0;

            // Compute the sum of weights. Need to do this first to have ready for sensitivities
            // std::cout << "[ ";    // brendan delete prints
            for ( uint iR = 0; iR < tNumRays; iR++ )
            {
                if ( std::abs( tRaysOnFacet( iR ).second - tMedian ) <= tStdDeviation )
                {
                    tSumOfValidWeights += tConeDistances.mRayCones.mDirectionWeights( iR );
                    // std::cout << iR << ", "; // brendan delete prints
                }
            }
            // std::cout << "], "; // brendan delete prints

            MORIS_ERROR( tSumOfValidWeights > 0.0, "Sum of valid weights is zero for facet %d. Check ray directions and intersections.", iF );

            // Get the vertex indices for this facet, these vertices need their sensitivities updated
            Vector< moris_index > tFacetVertices = this->get_facets_vertex_indices( iF );

            // Loop over ray intersections for this facet
            for ( uint iR = 0; iR < tNumRays; iR++ )
            {
                // Check if the ray distance is within one standard deviation from the median
                if ( std::abs( tRaysOnFacet( iR ).second - tMedian ) <= tStdDeviation )
                {
                    // Add the weighted distance to the weighted sum
                    tDiamWeightedSum += tRaysOnFacet( iR ).second * tConeDistances.mRayCones.mDirectionWeights( iR );

                    // Each ray intersection has sensitivities wrt to 2*dim vertices. All vertices of the origin facet and all vertices of the intersected facet
                    // Note that if a ray hits the neighboring facet, the origin facet vertices and the intersected facet vertices can share some vertices
                    // We compute the derivative of the shape diameter to the vertices and store the accumulation of sensitivity for all vertices
                    // For each vertex ddiameter/dv = dr/dx * weight / sum_of_weights

                    Vector< moris_index > tIntersectedFacetVertices = this->get_facets_vertex_indices( tRaysOnFacet( iR ).first );

                    // Origin facets vertices
                    for ( uint iV : tFacetVertices )
                    {
                        Matrix< DDRMat > tRotation = -trans( sdf::rotation_matrix( tConeDistances.mRayCones.mTheta( iR ) ) );

                        // Raycast distance is a function of the origin and the normal, which are both influenced by these vertices
                        Matrix< DDRMat > tCentroidSens = this->compute_draycast_dorigin( tCentroids.get_column( iF ), tConeDistances.mRayCones.mRayDirections( iF ).get_column( iR ), tRaysOnFacet( iR ).first )
                                                       * this->compute_dfacet_centroid_dvertex( iF, iV, true );
                        Matrix< DDRMat > tNormalSens = this->compute_draycast_ddirection( tCentroids.get_column( iF ), tConeDistances.mRayCones.mRayDirections( iF ).get_column( iR ), tRaysOnFacet( iR ).first )
                                                     * tRotation * this->compute_dfacet_normal_dvertex( iF, iV, true );

                        // Full sensitivity for this vertex
                        Matrix< DDRMat > tVertexSensitivity = tFacetMeasures( iF ) * tConeDistances.mRayCones.mDirectionWeights( iR ) * ( tCentroidSens + tNormalSens ) / tSumOfValidWeights;

                        // Accumulate the sensitivity for this vertex
                        mdShapeDiameterdVertex.set_row( iV, mdShapeDiameterdVertex.get_row( iV ) + tVertexSensitivity );
                    }

                    Matrix< DDRMat > tVertexSens = this->compute_draycast_dvertices( tCentroids.get_column( iF ), tConeDistances.mRayCones.mRayDirections( iF ).get_column( iR ), tRaysOnFacet( iR ).first );
                    for ( uint iV = 0; iV < tIntersectedFacetVertices.size(); iV++ )
                    {
                        // Raycast distance is a function of the facets vertices explicitly for these vertices. Computed for all vertices at once above
                        Matrix< DDRMat > tVertexSensitivity = tFacetMeasures( iF ) * tConeDistances.mRayCones.mDirectionWeights( iR ) * tVertexSens.get_row( iV ) / tSumOfValidWeights;

                        mdShapeDiameterdVertex.set_row( tIntersectedFacetVertices( iV ), mdShapeDiameterdVertex.get_row( tIntersectedFacetVertices( iV ) ) + tVertexSensitivity );
                    }
                }
            }

            // Divide by the sum of weights to get the weighted average
            tShapeDiameterValues( iF ) = tDiamWeightedSum / tSumOfValidWeights;
        }
        // std::cout << "];\n";    // brendan delete prints

        return tShapeDiameterValues;
    }

    // --------------------------------------------------------------------------------------------------------------

    real Surface_Mesh::compute_global_shape_diameter(
            real                   aConeAngle,
            uint                   aNumPolarRays,
            uint                   aNumAzimuthRays,
            Agglomeration_Function aAgglomerationFunction )
    {
        // Compute the shape diameter for every vertex in the surface mesh
        mShapeDiameters = this->compute_nodal_shape_diameter( aConeAngle, aNumPolarRays, aNumAzimuthRays );

        // Compute the area of every facet in the surface mesh
        Vector< real > tFacetMeasure     = this->compute_facet_measure();
        real           tTotalSurfaceArea = std::accumulate( tFacetMeasure.begin(), tFacetMeasure.end(), 0.0 );

        // Need to integrate the shape diameter over the surface mesh, loop over facets
        real tIntegratedShapeDiameter = 0.0;
        for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
        {
            // Add the shape diameter at this facet times the area of this facet to the integrated shape diameter
            tIntegratedShapeDiameter += mShapeDiameters( iF ) * tFacetMeasure( iF );
        }

        // Divide by the total area of the surface mesh to get the integrated shape diameter
        return mGlobalShapeDiameter = tIntegratedShapeDiameter / tTotalSurfaceArea;
    }

    // --------------------------------------------------------------------------------------------------------------

    real Surface_Mesh::compute_facet_measure( const uint aFacetIndex ) const
    {
        // Get the vertex coordinates of the facet
        Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( aFacetIndex );

        switch ( this->get_spatial_dimension() )
        {
            case 2:
            {
                // Length is the norm of the edge
                return norm( tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 ) );
                break;
            }
            case 3:
            {
                // Area is the norm of the cross product of two edges divided by 2
                Matrix< DDRMat > tEdge1 = tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 );
                Matrix< DDRMat > tEdge2 = tVertexCoordinates.get_column( 2 ) - tVertexCoordinates.get_column( 0 );
                return 0.5 * norm( cross( tEdge1, tEdge2 ) );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_facet_measure - Only implemented for 2D and 3D." );
                return 0.0;
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Vector< real > Surface_Mesh::compute_facet_measure() const
    {
        uint tNumCells = this->get_number_of_facets();

        Vector< real > tFacetMeasure( tNumCells );

        switch ( this->get_spatial_dimension() )
        {
            case 2:
            {
                for ( uint iF = 0; iF < tNumCells; iF++ )
                {
                    Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( iF );
                    tFacetMeasure( iF )                 = norm( tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 ) );
                }
                break;
            }
            case 3:
            {
                for ( uint iF = 0; iF < tNumCells; iF++ )
                {
                    Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( iF );

                    Matrix< DDRMat > tEdge1 = tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 );
                    Matrix< DDRMat > tEdge2 = tVertexCoordinates.get_column( 2 ) - tVertexCoordinates.get_column( 0 );
                    tFacetMeasure( iF )     = 0.5 * norm( cross( tEdge1, tEdge2 ) );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_facet_measure - Only implemented for 2D and 3D." );
                return tFacetMeasure;
            }
        }


        return tFacetMeasure;
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_dfacet_measure_dvertex( uint aFacetVertex, uint aVertexIndex, bool aRequireIsMember ) const
    {
        uint tDim = this->get_spatial_dimension();

        switch ( tDim )
        {
            case 2:
            {
                // Get the vertex indices of the facet
                const Vector< moris_index >& tFacetVertices = this->get_facets_vertex_indices( aFacetVertex );

                // Determine which index of the facet the requested vertex is, set sign accordingly
                real tSign = aVertexIndex == (uint)tFacetVertices( 0 ) ? 1.0 : ( aVertexIndex == (uint)tFacetVertices( 1 ) ? -1.0 : 0.0 );

                MORIS_ERROR( not aRequireIsMember or tSign != 0.0,
                        "Surface_Mesh::compute_dfacet_measure_dvertex - Vertex %d is not a member of facet %d and the function was called requiring this to be the case.",
                        aVertexIndex,
                        aFacetVertex );

                // Edge vector divided by the measure (signed based on the local vertex index within the facet)
                return tSign * trans( this->get_vertex_coordinates( tFacetVertices( 0 ) ) - this->get_vertex_coordinates( tFacetVertices( 1 ) ) ) / this->compute_facet_measure( aFacetVertex );
                break;
            }
            case 3:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_dfacet_measure_dvertex - Not implemented for 3D yet." );
                return { { 0.0, 0.0, 0.0 } };    // TODO: implement for 3D
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_dfacet_measure_dvertex - Only implemented for 2D and 3D." );
                return { {} };
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_ddiameter_dvertex(
            const uint aVertexIndex ) const
    {
        // FIXME BRENDAN add agglomeration function and sensitivity here

        // Compute the facet measures
        Vector< real > tFacetMeasures = this->compute_facet_measure();
        real           tSurfaceArea   = std::accumulate( tFacetMeasures.begin(), tFacetMeasures.end(), 0.0 );

        // Compute sum(l_i * d_i)
        real tWeightedDiameterSum = 0.0;
        for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
        {
            tWeightedDiameterSum += tFacetMeasures( iF ) * mShapeDiameters( iF );
        }

        // Get the facets connected to this vertex
        const Vector< moris_index >& tConnectedFacets = this->get_vertexs_facet_indices( aVertexIndex );

        Matrix< DDRMat > tIntegratedSensitivity( 1, this->get_spatial_dimension(), 0.0 );

        // Loop over connected facets to compute their contribution to the sensitivity as their area is changing
        for ( moris_index iF : tConnectedFacets )
        {
            Matrix< DDRMat > tdMeasuredVertex = this->compute_dfacet_measure_dvertex( iF, aVertexIndex, true );
            tIntegratedSensitivity += tdMeasuredVertex * mShapeDiameters( iF ) / tSurfaceArea - tdMeasuredVertex * tWeightedDiameterSum / ( tSurfaceArea * tSurfaceArea );
        }

        // Add shape diameter sensitivity for this vertex (precomputed during shape diameter computation)
        tIntegratedSensitivity += mdShapeDiameterdVertex.get_row( aVertexIndex ) / tSurfaceArea;

        return tIntegratedSensitivity;
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_dvolume_dvertex( uint aVertexIndex ) const
    {
        uint tDim = this->get_spatial_dimension();

        Matrix< DDRMat > tDVolumeDVertex( 1, tDim, MORIS_REAL_MAX );

        if ( this->get_spatial_dimension() == 2 )
        {
            // Get the vertices connected to this vertex
            const Vector< moris_index >& tNeighbors = mVertexToVertexConnectivity( aVertexIndex );

            moris_index tLowNeighbor;
            moris_index tHighNeighbor;
            if ( aVertexIndex > 1 and std::any_of( tNeighbors.begin(), tNeighbors.end(), []( moris_index i ) { return i == 0; } ) )
            {
                tLowNeighbor  = tNeighbors( 0 ) == 0 ? tNeighbors( 1 ) : tNeighbors( 0 );
                tHighNeighbor = 0;
            }
            else
            {    // Determine the order of the vertices
                tLowNeighbor  = tNeighbors( 0 ) < tNeighbors( 1 ) ? tNeighbors( 0 ) : tNeighbors( 1 );
                tHighNeighbor = tNeighbors( 0 ) > tNeighbors( 1 ) ? tNeighbors( 0 ) : tNeighbors( 1 );
            }

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

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_dfacet_centroid_dvertex( uint aFacet, uint aVertex, bool aRequireIsMember ) const
    {
        uint tDim = this->get_spatial_dimension();

        // Get the vertex indices of the facet
        const Vector< moris_index >& tFacetVertices = this->get_facets_vertex_indices( aFacet );

        // Get the index of the vertex in the facet
        uint tVertexInFacetIndex = MORIS_UINT_MAX;
        if ( aVertex == (uint)tFacetVertices( 0 ) )
        {
            tVertexInFacetIndex = 0;
        }
        else if ( aVertex == (uint)tFacetVertices( 1 ) )
        {
            tVertexInFacetIndex = 1;
        }
        else if ( tDim == 3 && aVertex == (uint)tFacetVertices( 2 ) )
        {
            tVertexInFacetIndex = 2;
        }

        // If the vertex is not part of the facet, return or throw an error based on the desired behavior
        if ( tVertexInFacetIndex == MORIS_UINT_MAX )
        {
            MORIS_ERROR( not aRequireIsMember, "Surface_Mesh::compute_dfacet_normal_dvertex - Vertex is not part of the facet and the function was called requiring this to be the case." );
            return Matrix< DDRMat >( tDim, tDim, 0.0 );
        }
        else
        {
            // Each vertex contributes equally to the centroid
            return ( 1.0 / static_cast< real >( tDim ) ) * eye( tDim, tDim );
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_dfacet_normal_dvertex( uint aFacet, uint aVertex, bool aRequireIsMember ) const
    {
        uint tDim = this->get_spatial_dimension();

        // Initialize sensitivity
        Matrix< DDRMat > tSensitivity( tDim, tDim, 0.0 );

        // Get the vertex indices of the facet
        const Vector< moris_index >& tFacetVertices = this->get_facets_vertex_indices( aFacet );

        // Get the index of the vertex in the facet
        uint tVertexInFacetIndex = MORIS_UINT_MAX;
        if ( aVertex == (uint)tFacetVertices( 0 ) )
        {
            tVertexInFacetIndex = 0;
        }
        else if ( aVertex == (uint)tFacetVertices( 1 ) )
        {
            tVertexInFacetIndex = 1;
        }
        else if ( tDim == 3 && aVertex == (uint)tFacetVertices( 2 ) )
        {
            tVertexInFacetIndex = 2;
        }

        // If the vertex is not part of the facet, return or throw an error based on the desired behavior
        if ( tVertexInFacetIndex == MORIS_UINT_MAX )
        {
            MORIS_ERROR( not aRequireIsMember, "Surface_Mesh::compute_dfacet_normal_dvertex - Vertex is not part of the facet and the function was called requiring this to be the case." );
            return Matrix< DDRMat >( tDim, tDim, 0.0 );
        }

        // Get the vertex coordinates of the facet
        Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( aFacet );

        switch ( tDim )
        {
            case 2:
            {
                // magnitude of the normal vector
                real tInverseNormalVectorNorm = 1.0 / norm( tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 ) );

                if ( tVertexInFacetIndex == 0 )
                {
                    tSensitivity = { { ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 0, 0 ) ) * ( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 1, 0 ) ), -1.0 * std::pow( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 0, 0 ), 2.0 ) }, { std::pow( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 1, 0 ), 2.0 ), ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 0, 0 ) ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 1, 1 ) ) } };
                }
                else
                {
                    tSensitivity = { { ( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 1, 0 ) ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 0, 1 ) ), std::pow( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 0, 0 ), 2.0 ) }, { -1.0 * std::pow( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 1, 0 ), 2.0 ), ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 0, 1 ) ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 1, 1 ) ) } };
                }
                return std::pow( tInverseNormalVectorNorm, 3.0 ) * tSensitivity;
            }
            case 3:
            {
                // Compute the normal vector (not unit)
                Matrix< DDRMat > tNormal = cross( tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 ), tVertexCoordinates.get_column( 2 ) - tVertexCoordinates.get_column( 0 ) );

                // magnitude of the normal vector
                real tInverseNormalVectorNorm = 1.0 / norm( tNormal );

                Matrix< DDRMat > tNormalVectorNormSensitivity( 3, 3, 0.0 );

                if ( tVertexInFacetIndex == 0 )
                {
                    tSensitivity                 = { { 0.0, tVertexCoordinates( 2, 1 ) - tVertexCoordinates( 2, 2 ), tVertexCoordinates( 1, 2 ) - tVertexCoordinates( 1, 1 ) }, { tVertexCoordinates( 2, 2 ) - tVertexCoordinates( 2, 1 ), 0.0, tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 0, 2 ) }, { tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 1, 2 ), tVertexCoordinates( 0, 2 ) - tVertexCoordinates( 0, 1 ), 0.0 } };
                    tNormalVectorNormSensitivity = { { tNormal( 2 ) * ( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 1, 2 ) ) - tNormal( 1 ) * ( tVertexCoordinates( 2, 1 ) - tVertexCoordinates( 2, 2 ) ), -tNormal( 2 ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 0, 2 ) ) + tNormal( 0 ) * ( tVertexCoordinates( 2, 1 ) - tVertexCoordinates( 2, 2 ) ), tNormal( 1 ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 0, 2 ) ) - tNormal( 0 ) * ( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 1, 2 ) ) } };
                }
                else if ( tVertexInFacetIndex == 1 )
                {
                    tSensitivity                 = { { 0.0, tVertexCoordinates( 2, 2 ) - tVertexCoordinates( 2, 0 ), tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 1, 2 ) }, { tVertexCoordinates( 2, 0 ) - tVertexCoordinates( 2, 2 ), 0.0, tVertexCoordinates( 0, 2 ) - tVertexCoordinates( 0, 0 ) }, { tVertexCoordinates( 1, 2 ) - tVertexCoordinates( 1, 0 ), tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 0, 2 ), 0.0 } };
                    tNormalVectorNormSensitivity = { { -tNormal( 2 ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 1, 2 ) ) + tNormal( 1 ) * ( tVertexCoordinates( 2, 0 ) - tVertexCoordinates( 2, 2 ) ), tNormal( 2 ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 0, 2 ) ) - tNormal( 0 ) * ( tVertexCoordinates( 2, 0 ) - tVertexCoordinates( 2, 2 ) ), -tNormal( 1 ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 0, 2 ) ) + tNormal( 0 ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 1, 2 ) ) } };
                }
                else if ( tVertexInFacetIndex == 2 )
                {
                    tSensitivity                 = { { 0.0, tVertexCoordinates( 2, 0 ) - tVertexCoordinates( 2, 1 ), tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 1, 0 ) }, { tVertexCoordinates( 2, 1 ) - tVertexCoordinates( 2, 0 ), 0.0, tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 0, 1 ) }, { tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 1, 1 ), tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 0, 0 ), 0.0 } };
                    tNormalVectorNormSensitivity = { { tNormal( 2 ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 1, 1 ) ) - tNormal( 1 ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 2, 1 ) ), -tNormal( 2 ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 0, 1 ) ) + tNormal( 0 ) * ( tVertexCoordinates( 2, 0 ) - tVertexCoordinates( 2, 1 ) ), tNormal( 1 ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 0, 1 ) ) - tNormal( 0 ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 1, 1 ) ) } };
                }

                return tInverseNormalVectorNorm * tInverseNormalVectorNorm * ( tSensitivity - tNormal * tInverseNormalVectorNorm * tNormalVectorNormSensitivity );
            }
            default:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_dfacet_normal_dvertex - Only implemented for 2D and 3D." );
                return { {} };
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_draycast_dorigin( const Matrix< DDRMat >& aOrigin, const Matrix< DDRMat >& aDirection, uint aFacetIndex ) const
    {
        uint tDim = this->get_spatial_dimension();

        // Get the vertex coordinates of the facet
        Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( aFacetIndex );

        switch ( tDim )
        {
            case 2:
            {
                // Compute the edge vector
                Matrix< DDRMat > tEdge = tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 );

                // Compute the inverse determinant
                real tInvDet = 1.0 / cross_2d( aDirection, tEdge );

                MORIS_ASSERT( std::abs( tInvDet ) < MORIS_REAL_MAX, "Surface_Mesh::compute_draycast_dorigin - Division by zero detected in inverse determinant computation." );

                Matrix< DDRMat > tNormal = { { -tEdge( 1 ), tEdge( 0 ) } };    // non-unit vector. Transposed to make the sensitivity matrix the correct size

                return tInvDet * tNormal;
            }
            case 3:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_draycast_dorigin - do 3d implementation" );
                return { {} };
            }
            default:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_draycast_dorigin - Only implemented for 2D and 3D." );
                return { {} };
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_draycast_ddirection( const Matrix< DDRMat >& aOrigin, const Matrix< DDRMat >& aDirection, uint aFacetIndex ) const
    {
        uint tDim = this->get_spatial_dimension();

        // Get the vertex coordinates of the facet
        Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( aFacetIndex );

        switch ( tDim )
        {
            case 2:
            {
                // Compute the edge vector
                Matrix< DDRMat > tEdge = tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 );

                // Vector from the origin of the ray to the origin of the first vertex
                Matrix< DDRMat > tOriginToV0;
                if ( aOrigin.n_cols() == 1 )
                {
                    tOriginToV0 = tVertexCoordinates.get_column( 0 ) - aOrigin;
                }
                else
                {
                    tOriginToV0 = tVertexCoordinates.get_column( 0 ) - trans( aOrigin );
                }

                // Compute the determinant of the ray to the first vertex and the edge
                real tDet = cross_2d( tOriginToV0, tEdge );

                // Compute the inverse determinant of the edge and the ray direction
                real tInvDet = 1.0 / cross_2d( aDirection, tEdge );

                Matrix< DDRMat > tNormal = { { -tEdge( 1 ), tEdge( 0 ) } };    // non-unit vector. Transposed to make the sensitivity matrix the correct size

                return tInvDet * tInvDet * tDet * tNormal;
            }
            case 3:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_draycast_dorigin - do 3d implementation" );
                return { {} };
            }
            default:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_draycast_dorigin - Only implemented for 2D and 3D." );
                return { {} };
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_draycast_dvertices( const Matrix< DDRMat >& aOrigin, const Matrix< DDRMat >& aDirection, uint aFacetIndex ) const
    {
        uint tDim = this->get_spatial_dimension();

        Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates_of_facet( aFacetIndex );

        switch ( tDim )
        {
            case 2:
            {
                // Compute the edge vector
                Matrix< DDRMat > tEdge = tVertexCoordinates.get_column( 1 ) - tVertexCoordinates.get_column( 0 );

                // Vector from the origin of the ray to the origin of the first vertex
                Matrix< DDRMat > tOriginToV0;
                if ( aOrigin.n_cols() == 1 )
                {
                    tOriginToV0 = tVertexCoordinates.get_column( 0 ) - aOrigin;
                }
                else
                {
                    tOriginToV0 = tVertexCoordinates.get_column( 0 ) - trans( aOrigin );
                }

                // Compute the inverse determinant of the edge and the ray direction
                real tInvDet = 1.0 / cross_2d( aDirection, tEdge );

                // Determinant of direction and origin
                real tDirOriginDet = cross_2d( aOrigin, aDirection );

                // Determinant of direction and vertices
                real tDirV1Det = cross_2d( tVertexCoordinates.get_column( 1 ), aDirection );
                real tDirV0Det = cross_2d( tVertexCoordinates.get_column( 0 ), aDirection );

                Matrix< DDRMat > tNormal  = { { tEdge( 1 ), -tEdge( 0 ) }, { -tEdge( 1 ), tEdge( 0 ) } };                  // non-unit vector. Transposed to make the sensitivity matrix the correct size
                Matrix< DDRMat > tScaling = { { tDirOriginDet - tDirV1Det, 0.0 }, { 0.0, tDirOriginDet - tDirV0Det } };    // Determinant scaling per vertex

                return tScaling * tNormal * tInvDet * tInvDet;
            }
            case 3:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_draycast_dvertices - do 3d implementation" );
                return { {} };
            }
            default:
            {
                MORIS_ERROR( false, "Surface_Mesh::compute_draycast_dorigin - Only implemented for 2D and 3D." );
                return { {} };
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Surface_Mesh::get_nodal_shape_diameter_sensitivities() const
    {
        return mdShapeDiameterdVertex;
    }

    // --------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh::compute_vertex_normals() const
    {
        const uint             tNumVertices  = this->get_number_of_vertices();
        const uint             tDim          = this->get_spatial_dimension();
        const Matrix< DDRMat > tFacetNormals = this->get_all_facet_normals();
        const Vector< real >   tFacetMeasure = this->compute_facet_measure();
        Matrix< DDRMat >       tVertexNormals( tDim, tNumVertices );

        auto tNormal = Matrix< DDRMat >( tDim, 1 );
        for ( moris::size_t iV = 0; iV < tNumVertices; iV++ )
        {
            Vector< moris_index > tVertexCellNeighbors = this->get_vertexs_facet_indices( iV );
            auto const            tNumNeighbors        = static_cast< moris::size_t >( tVertexCellNeighbors.size() );
            tNormal.fill( 0.0 );    // reset the current normal to zero for each vertex normal calculation

            // compute the normal as the weighted average of the facet normals of the neighboring facets
            for ( moris::size_t iVN = 0; iVN < tNumNeighbors; iVN++ )
            {
                int const tFacetIndex = tVertexCellNeighbors( iVN );
                tNormal += tFacetNormals.get_column( tFacetIndex ) * tFacetMeasure( tFacetIndex );
            }
            tVertexNormals.set_column( iV, tNormal / norm( tNormal ) );
        }

        return tVertexNormals;
    }

}    // namespace moris::mtk
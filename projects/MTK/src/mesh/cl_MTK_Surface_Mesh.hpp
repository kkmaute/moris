/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Surface_Mesh.hpp
 *
 */
#pragma once

#include "cl_MTK_Enums.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"

#if MORIS_HAVE_ARBORX
#include <ArborX.hpp>
#include <ArborX_Box.hpp>

#include <Kokkos_Macros.hpp>
#include <Kokkos_View.hpp>
#include <decl/Kokkos_Declare_SERIAL.hpp>
namespace moris::mtk::arborx
{
    template< typename MemorySpace >
    struct QueryRays;

    /**
     * @brief Converts a moris::Matrix< moris::DDRMat > to an ArborX::Point or ArborX::Vector.
     * @tparam T The type of the ArborX object to be returned (Point or Vector)
     * @param aMatrix The matrix to be converted (either 3x1 or 2x1)
     * @return The converted ArborX object
     */
    template< typename T >
    T coordinate_to_arborx_point( moris::Matrix< moris::DDRMat > const & aMatrix );
}    // namespace moris::mtk::arborx

namespace moris::mtk
{
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using MemorySpace    = ExecutionSpace::memory_space;
}    // namespace moris::mtk


#endif

namespace moris::mtk
{

    // helper function for 2d raycast
    real cross_2d( const Matrix< DDRMat >& aVector1, const Matrix< DDRMat >& aVector2 );

    typedef Vector< std::pair< uint, real > > Intersection_Vector;    // pair of (facet index, distance)

    struct Ray_Cones
    {
        Vector< real > mDirectionWeights;    // Weight for each ray direction

        Vector< Matrix< DDRMat > > mRayDirections;    // Each entry in the vector corresponds to a vertex. Each Matrix contains the ray directions as columns.

        /**
         * Constructor to size the data structures
         */
        Ray_Cones( uint aNumVertices, uint aNumRays, uint aDim )
                : mDirectionWeights( aNumRays )
                , mRayDirections( aNumVertices, Matrix< DDRMat >( aDim, aNumRays ) )
        {
        }
    };

    struct Shape_Diameter_Distances
    {
        Vector< real >& mDirectionWeights;    // Weight for each ray direction

        // Each entry in the vector corresponds to a vertex.
        // Each vertex has a number of rays( aNumPolarRays* aNumAzimuthRays ) associated with it.
        // Each ray only stores its nearest intersection (facet index, distance)
        // Usage: Input( vertex index, ray index ) -> Output( facet index, distance )
        Vector< Intersection_Vector > mDistances;

        /**
         * Constructor to size the data structures
         */
        Shape_Diameter_Distances( Vector< real >& aDirectionWeights, uint aNumVertices, uint aNumRays )
                : mDirectionWeights( aDirectionWeights )
                , mDistances( aNumVertices, Intersection_Vector( aNumRays ) )
        {
        }
    };

    /**
     * @brief This class is used to extract a surface mesh from a given mesh.
     * This class does not store any vertex or cell data but only provides the necessary information to access the data in the mesh.
     * The base surface mesh will ONLY use local indices that only refer to the vertices/facets in the surface mesh.
     */
    class Surface_Mesh
    {
      public:
        /**
         * @brief Constructor for base surface mesh class. Builds mVertexToCellIndices, mFacetNormals, and mVertexCoordinates.
         *
         */
        Surface_Mesh( const Matrix< DDRMat >&          aVertexCoordinates,
                const Vector< Vector< moris_index > >& aFacetConnectivity,
                real                                   aIntersectionTolerance = 1e-8 );

        /**
         * Removes vertices that are not part of any facet and updates the facet connectivity accordingly
         *
         * @param aVertexCoordinates <dimension> x <number of vertices> matrix containing vertex coordinates with potential extraneous vertices
         * @param aFacetConnectivity Vector of vectors containing the local indices of the vertices that form each facet
         */
        void clean_extraneous_vertices();

        // -------------------------------------------------------------------------------
        // Mesh deformation methods
        // -------------------------------------------------------------------------------

        /**
         * Sets the displacements for all facet vertices at once and updates the normal vectors accordingly
         *
         * @param aDisplacements <dimension> x <number of vertices> matrix containing displacment data for all vertices
         */
        virtual void set_all_displacements( const Matrix< DDRMat >& aDisplacements );

        /**
         * Sets the displacement for ONE vertex
         * NOTE: This method does NOT update the facet normals. Use set_all_displacements() to update the normals or ensure that the normals are updated manually by calling initialize_facet_normals()
         *
         * @param aVertexIndex vertex index to set displacement for
         * @param aDisplacement <dimension > x < 1 > matrix containing vertex's displacment
         */
        void set_vertex_displacement( const uint aVertexIndex, const Matrix< DDRMat >& aDisplacement );

        /**
         * @brief removes any displacement, rotation, and scaling applied to the surface mesh
         */
        void reset_coordinates();

        // -------------------------------------------------------------------------------
        // Accessor methods
        // -------------------------------------------------------------------------------

        [[nodiscard]] virtual const Matrix< DDRMat > get_all_vertex_coordinates() const;

        /**
         * @brief Gets the coordinates of a single vertex from the local index aVertexIndex
         */
        [[nodiscard]] virtual const Matrix< DDRMat > get_vertex_coordinates( const uint aVertexIndex ) const;

        /**
         * @brief Gets the original coordinates (no displacement added) of all vertices in the surface mesh
         * Size: < spatial dim x number of vertices >
         */
        [[nodiscard]] virtual const Matrix< DDRMat > get_all_original_vertex_coordinates() const;

        /**
         * Gets the original coordinates of a single vertex from the local index aVertexIndex
         */
        [[nodiscard]] virtual const Matrix< DDRMat > get_original_vertex_coordinates( const uint aVertexIndex ) const;

        /**
         * @brief gets the displacements of all vertices in the surface mesh
         */
        [[nodiscard]] virtual const Matrix< DDRMat >& get_vertex_displacements() const;

        /**
         * @brief gets the entire facet to vertex connectivity of the surface mesh
         */
        [[nodiscard]] const Vector< Vector< moris_index > >& get_facet_connectivity() const;

        /**
         * @brief gets the entire vertex to facet connectivity of the surface mesh
         */
        [[nodiscard]] const Vector< Vector< moris_index > >& get_vertex_connectivity() const;

        /**
         * @brief Gets the indices to the vertices that form the facet with the local index aFacetIndex
         *
         * @param aFacetIndex local index of the facet
         * @return Vector< moris_index > local vertex indices that form the facet
         */
        [[nodiscard]] const Vector< moris_index >& get_facets_vertex_indices( const uint aFacetIndex ) const;

        /**
         * @brief Gets the indices to the facets that are connected to the vertex with the local index aVertexIndex
         *
         * @param aVertexIndex local index of the vertex
         * @return Vector< moris_index > local facet indices that are connected to the vertex
         */
        [[nodiscard]] const Vector< moris_index >& get_vertexs_facet_indices( const uint aVertexIndex ) const;

        /**
         * @brief Gets the index of all vertices connected by a facet to the vertex with the local index aVertexIndex
         *
         * @param aVertexIndex local index of the vertex
         * @return Vector< moris_index > local vertex indices that are connected to the vertex through a facet
         */
        [[nodiscard]] const Vector< moris_index >& get_vertex_neighbors( const uint aVertexIndex ) const;

        /**
         * @brief Gets the coordinates of all vertices that form the facet with the local index aFacetIndex
         * Size: < spatial dim x number of vertices in the facet >
         */
        [[nodiscard]] Matrix< DDRMat > get_all_vertex_coordinates_of_facet( const uint aFacetIndex ) const;

        /**
         * @brief Returns the facet normals for each facet in the surface mesh.
         * @return A (d x n) matrix where d is the dimension of the mesh and n is the number of facets in the surface mesh (holding the normal components).
         */
        [[nodiscard]] const Matrix< DDRMat >& get_all_facet_normals() const;

        /**
         * @brief Gets the normal vector of the facet with the local index aFacetIndex
         *
         * @param aFacetIndex local index of the facet
         * @return Matrix< DDRMat > normal vector of the facet
         */
        [[nodiscard]] const Matrix< DDRMat > get_facet_normal( const uint aFacetIndex ) const;

        [[nodiscard]] virtual uint get_spatial_dimension() const;

        /**
         * @brief Returns the number of facets in the mesh
         */
        [[nodiscard]] uint get_number_of_facets() const;

        [[nodiscard]] uint get_number_of_vertices() const;

        [[nodiscard]] real get_intersection_tolerance() const;


        // -------------------------------------------------------------------------------
        // Raycast methods
        // -------------------------------------------------------------------------------

        /**
         * @brief Determines if a point is inside or outside the surface mesh via raycasting
         * This method utilizes ArborX to find ray facet intersections, and then computes the intersection locations for the ray.
         * The region is determined by the number of intersections. Even number = outside, Odd number = inside.
         *
         * @param aPoint Ray origin point
         */
        Mesh_Region
        get_region_from_raycast( const Matrix< DDRMat >& aPoint ) const;

        /**
         * @brief Determines if the points are inside or outside the surface mesh via raycasting
         * This method utilizes ArborX to find ray facet intersections, and then computes the intersection locations for the ray.
         * The region is determined by the number of intersections. Even number = outside, Odd number = inside.
         *
         * @param aPoint Ray origin points. Each column is a point, size <dimension> x <number of points>.
         */
        Vector< Mesh_Region >
        batch_get_region_from_raycast( Matrix< DDRMat >& aPoint ) const;

        /**
         * @brief Determines if a point is inside or outside the surface mesh via raycasting
         * This method utilizes ArborX to find ray facet intersections, and then computes the intersection locations for the ray.
         * The method will only cast a single ray, which may not be sufficient to determine the region if the ray hits an edge or another pathological case is detected
         * The region is determined by the number of intersections. Even number = outside, Odd number = inside.
         *
         * @param aPoint Ray origin point. Passed by value as it may be altered
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @param aWarning Flag if the ray hits an edge or another pathological case
         * @param aIgnoreWarnings If true, all rays that hit a warning will be ceased immediately. Otherwise, the ray is cast and the result is returned as normal.
         */
        Intersection_Vector
        cast_single_ray(
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection,
                bool&                   aWarning,
                bool                    aIgnoreWarnings = true ) const;

        /**
         * Casts many rays and returns all of the associated intersection pairs. Uses the same directions for every origin.
         *
         * @param aOrigins Origin points for the rays. Each column is a new origin, size <dimension> x <number of origins>
         * @param aDirections Directions for the rays, where each column is a new direction. Size <dimension> x <number of directions>
         * @param aWarnings Vector of warnings for each ray. Size <number of origins>
         * @param aIgnoreWarnings If true, all rays that hit a warning will be ceased immediately. Otherwise, the ray is cast and the result is returned as normal.
         * @return Vector< Vector< Intersection_Vector > > Outer vector corresponds to the origin point, and the inner vector corresponds to the direction.
         *      Intersection_Vector contains pairs of <uint, real> which correspond to the facet index and the distance to the intersection point
         */
        Vector< Vector< Intersection_Vector > >
        cast_batch_of_rays(
                Matrix< DDRMat >&         aOrigins,
                Matrix< DDRMat >&         aDirections,
                Vector< Vector< bool > >& aWarnings,
                bool                      aIgnoreWarnings = true ) const;

        /**
         * Casts many rays and returns all of the associated intersection pairs. Allows for raycasting for any direction for any of the points.
         *
         * @param aOrigins Origin points for the rays. Each column is a new origin, size <dimension> x <number of origins>
         * @param aDirections Directions for the rays. The size of the vector is <number of origins >,
         *      and each matrix in the vector contains directions, where each column is a new direction. Size <dimension> x <number of directions> )
         * @param aWarnings Vector of warnings for each ray. Size <number of origins>
         * @param aIgnoreWarnings If true, all rays that hit a warning will be ceased immediately. Otherwise, the ray is cast and the result is returned as normal.
         * @return Vector< Vector< Intersection_Vector > > Outer vector corresponds to the origin point, and the inner vector corresponds to the direction.
         *      Intersection_Vector contains pairs of <uint, real> which correspond to the facet index and the distance to the intersection point
         */
        Vector< Vector< Intersection_Vector > >
        cast_batch_of_rays(
                Matrix< DDRMat >&           aOrigins,
                Vector< Matrix< DDRMat > >& aDirections,
                Vector< Vector< bool > >&   aWarnings,
                bool                        aIgnoreWarnings = true ) const;

        //-------------------------------------------------------------------------------
        // Quantities of interest
        // -------------------------------------------------------------------------------

        /**
         * Computes the volume enclosed by the surface mesh
         */
        real compute_volume() const;

        /**
         * Gets the derivative of the volume wrt a vertex's coordinates
         *
         * @param aVertexIndex local index of the vertex to get sensitivities of
         * @return Matrix< DDRMat > dVolume/dVertex. Size: <1> x <spatial dim>
         */
        Matrix< DDRMat > compute_dvolume_dvertex( const uint aVertexIndex ) const;

        /**
         * Constructs a cone of rays for each vertex, centered around the vertex normal
         * @param aConeAngle Angle of the cone in degrees
         * @param aNumPolarRays Number of rays in the polar direction (θ)
         * @param aNumAzimuthRays Number of rays in the azimuth direction (φ) 1 if 2D
         * @return Ray_Cones struct which holds the weights and the ray directions for each vertex
         */
        Ray_Cones build_ray_cone_directions( real aConeAngle, uint aNumPolarRays, uint aNumAzimuthRays = 1 ) const;

        /**
         * Casts a cone of rays from each node, centered around the vertex normal, and gets the distances for each ray-facet intersection
         *
         * @param aConeAngle Angle of the cone in degrees
         * @param aNumPolarRays Number of rays in the polar direction (θ)
         * @param aNumAzimuthRays Number of rays in the azimuth direction (φ) 1 if 2D
         * @return Shape_Diameter_Distances struct which holds the weights and the intersection distances for each ray at each vertex
         */
        Shape_Diameter_Distances compute_nodal_shape_diameter_distances( real aConeAngle, uint aNumPolarRays, uint aNumAzimuthRays = 1 ) const;

        /**
         * @brief For each raycast result, determines the nearest intersection that is not part of the originating vertex's facets
         * Used for shape diameter computation
         *
         * @param aAllDistances Raycast results for ray for each vertex - result of a batch raycast
         * @return Vector< Intersection_Vector > Nearest non-trivial intersection for each ray at each vertex. Size: <number of vertices> x <number of rays per vertex>
         */
        Vector< Intersection_Vector > determine_nearest_nontrivial_intersections( Vector< Vector< Intersection_Vector > >& aAllIntersections ) const;

        /**
         * Computes the shape diameter for each surface mesh node
         * The shape diameter is computed by casting a cone of rays from each node and taking the minimum ray length of all the rays
         */
        Vector< real > compute_nodal_shape_diameter();

        /**
         * Computes the norm(brendan ?) of the shape diameter for every surface mesh node
         * The shape diameter is computed by casting a cone of rays from each node and taking the minimum ray length of all the rays
         *
         * @return real Norm of the shape diameter for each node
         */
        real compute_shape_diameter();

        /**
         * Computes the centroids of all the facets in the surface mesh
         * Size: <spatial dim> x <number of facets>
         * 
         * @return Matrix< DDRMat > facet centroids (average of the vertex coordinates of each facet)
         */
        Matrix< DDRMat > compute_facet_centroids() const;
        
        /**
         * Gets the derivative of the shape diameter wrt a vertex's coordinates
         *
         * @param aVertexIndex local index of the vertex to get sensitivities of
         * @return Matrix< DDRMat > dDiameter/dVertex. Size: <1> x <spatial dim>
         */
        Matrix< DDRMat > compute_ddiameter_dvertex( const uint aVertexIndex ) const;

        /**
         * Computes the measure (area in 3D, length in 2D) of a facet in the surface mesh
         * Virtual as some child implementations may prefer to compute and store this value. This implementation computes it on the fly.
         *
         * @param aFacetIndex local index of the facet to compute measure of.
         * @return Vector< real > facet measures. Size: <number of facets> x <1>
         */
        real compute_facet_measure( uint aFacetIndex ) const;

        /**
         * Computes the measure (area in 3D, length in 2D) of all facets in the surface mesh
         * Virtual as some child implementations may prefer to compute and store this value. This implementation computes it on the fly.
         *
         * @return Vector< real > facet measures. Size: <number of facets> x <1>
         */
        virtual Vector< real > compute_facet_measure() const;

        /**
         * @brief Gets the derivative of the facet measure wrt a vertex's coordinates
         *
         * @param aFacetIndex local index of the facet to get sensitivities of
         * @param aVertexIndex local index of the vertex to get sensitivities of
         * @return Matrix< DDRMat > dMeasure/dVertex. Size: <1> x <spatial dim>
         */
        virtual Matrix< DDRMat > compute_dfacet_measure_dvertex( const uint aFacetIndex, const uint aVertexIndex ) const;

        /**
         * Computes the normal vector for each vertex as the average of the normals of its facets
         * Virtual as some child implementations may prefer to compute and store this value. This implementation computes it on the fly.
         *
         * @return Matrix< DDRMat > vertex normals stored column-wise. Size: <spatial dim> x <number of vertices>
         */
        virtual Matrix< DDRMat > compute_vertex_normals() const;

        //-------------------------------------------------------------------------------
        // Output Methods
        // -------------------------------------------------------------------------------

        void write_to_file( const std::string& aFilePath ) const;

        //-------------------------------------------------------------------------------
        // Mesh modification methods
        // -------------------------------------------------------------------------------

      protected:
        /**
         * @brief Set the specified vertex's coordinates
         *
         * @param aVertexIndex local index of the vertex
         * @param aCoordinates desired coordinates for the vertex
         */
        void set_vertex_coordinates( const uint aVertexIndex, const Matrix< DDRMat >& aCoordinates );

        /**
         * @brief computes the facet normals and stores in mFacetNormals. mVertexCoordinates must be initialized first
         */
        void initialize_facet_normals();

        /**
         * @brief Determines which vertices are connected to each other via a facet. Does not store any information about the facets
         */
        void build_vertex_connectivity();

        /**
         * Determines which facets are connected to each vertex. Used for computing vertex normals
         */
        void build_vertex_to_facet_connectivity();

        //-------------------------------------------------------------------------------
        // Private methods useful for raycasting
        // -------------------------------------------------------------------------------

      private:
        /**
         * Computes the intersection location of a ray with a given facet.
         *
         * @param aFacet local facet index to compute intersection location
         * @param aPoint origin point of the ray
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @return real magnitude of the distance between the origin point and the intersection point. NaN if the facet is not intersected
         */
        real
        moller_trumbore(
                uint                    aFacet,
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection,
                bool&                   aWarning ) const;

        /**
         * @brief Moller trumbore algorithm for determining if the ray intersects a triangle.
         *
         * @param aFacet local facet index to compute intersection location
         * @param aAxis coordinate axis the ray is cast down
         * @param aPoint origin point of the ray
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @return real magnitude of the distance between the origin point and the intersection point. NaN if the facet is not intersected
         */
        real
        moller_trumbore_3D(
                uint                    aFacet,
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection,
                bool&                   aWarning ) const;

        /**
         * @brief (Psuedo) Moller trumbore algorithm for determining if the ray intersects a line.
         *
         * @param aFacet local facet index to compute intersection location
         * @param aPoint origin point of the ray
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @return real magnitude of the distance btween the origin point and the intersection point. NaN if the facet is not intersected
         */
        real
        moller_trumbore_2D(
                uint                    aFacet,
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection,
                bool&                   aWarning ) const;

        /**
         * @brief Takes candidate facets and attempts to compute the intersection locations for the given ray.
         *
         * @param aPoint Origin point of the ray
         * @param aDirection Direction of the ray
         * @param aCandidateFacets Indices of facets that the ray could hit.
         * @return Intersection locations and associated facet indices for every ray-facet intersection. Max size of aCandidateFacets.
         */
        Intersection_Vector determine_valid_intersections_from_candidates(
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection,
                const Vector< uint >&   aCandidateFacets,
                bool&                   aWarning,
                bool                    aIgnoreWarnings = true ) const;

        /**
         * @brief Removes duplicate intersections and sorts them by distance from the origin point
         *
         * @param aIntersections Intersection distances
         * @param aIntersectionFacetIndices Index of the facet that the ray hit.
         * @return Intersection_Vector Index, distance pairs sorted by distance
         */
        Intersection_Vector sort_and_find_unique_intersections( Intersection_Vector& aIntersections ) const;

        // Generates a random direction vector for raycasting
        Matrix< DDRMat > random_direction() const;

        //-------------------------------------------------------------------------------
        // ArborX API methods - used for raycasting acceleration
        // -------------------------------------------------------------------------------
#if MORIS_HAVE_ARBORX

        /**
         * Constructs the ArborX rays for the given points and directions
         *
         * @tparam MemorySpace
         * @tparam ExecutionSpace
         * @param aExecutionSpace
         * @param aOrigins Origin points for the rays. Each column is a new origin, size <dimension> x <number of origins>
         * @param aDirections Directions for each ray. Outer index corresponds to an origin, and each column of the matrix is a new direction.
         * Each origin can have its own number of directions. Size <number of origins > ( <dimension> x <number of directions> )
         * @return QueryRays< MemorySpace > Struct for ArborX ray queries
         */
        template< typename MemorySpace, typename ExecutionSpace >
        static arborx::QueryRays< MemorySpace > build_arborx_ray_batch(
                ExecutionSpace const &      aExecutionSpace,
                Matrix< DDRMat >&           aOrigins,
                Vector< Matrix< DDRMat > >& aDirections );

        /**
         * Constructs the ArborX rays for the given points and directions
         *
         * @tparam MemorySpace
         * @tparam ExecutionSpace
         * @param aExecutionSpace
         * @param aOrigins Origin points for the rays. Each column is a new origin, size <dimension> x <number of origins>
         * @param aDirections Directions for each ray. For this case, each origin will have the same directions. size <dimension> x <number of directions> )
         * @return QueryRays< MemorySpace > Struct for ArborX ray queries
         */
        template< typename MemorySpace, typename ExecutionSpace >
        static arborx::QueryRays< MemorySpace > build_arborx_ray_batch(
                ExecutionSpace const & aExecutionSpace,
                Matrix< DDRMat >&      aOrigins,
                Matrix< DDRMat >&      aDirections );

        /**
         * @brief Uses ArborX bounding volume hierarchy to determine which facets may be intersected by the ray
         * The ray only travels in the positive direction.
         *
         * @param aPoint Ray origin point
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @return Vector< uint > Local index of facets for which to check intersection
         */
        Vector< uint > preselect_with_arborx(
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection ) const;

        /**
         * @brief Uses ArborX bounding volume hierarchy to determine which facets may be intersected by the ray
         * The ray only travels in the positive direction.
         *
         * @param aOrigins Origin points for the rays. Each column is a new origin, size <dimension> x <number of origins>
         * @param aDirections Directions for each ray. Supports different directions for each origin point. size <number of origins > ( <dimension> x <number of directions> )
         * @return Vector< Vector< Vector< uint > > > Innermost vector is the indices of the facets that the ray hit. Middle vector is for each direction for each origin point. Outer vector is for origin points
         */
        Vector< Vector< Vector< uint > > > batch_preselect_with_arborx(
                Matrix< DDRMat >&           aOrigins,
                Vector< Matrix< DDRMat > >& aDirections ) const;

        /**
         * @brief Uses ArborX bounding volume hierarchy to determine which facets may be intersected by the ray
         * The ray only travels in the positive direction.
         *
         * @param aOrigins Origin points for the rays. Each column is a new origin, size <dimension> x <number of origins>
         * @param aDirections Directions for each ray. This version casts the same directions for every origin.size <dimension> x <number of directions>
         * @return Vector< Vector< Vector< uint > > > Innermost vector is the indices of the facets that the ray hit. Middle vector is for each direction for each origin point. Outer vector is for origin points
         */
        Vector< Vector< Vector< uint > > > batch_preselect_with_arborx(
                Matrix< DDRMat >& aOrigins,
                Matrix< DDRMat >& aDirections ) const;

      protected:
        /**
         * @brief Constructs the ArborX bounding volume hierarchy for the surface mesh and stores it in mBVH
         *
         */
        void construct_bvh();

#endif

        // -------------------------------------------------------------------------------
        // Member data
        // -------------------------------------------------------------------------------

      private:    // variables
        /**
         * @brief Stores the coordinates of all vertices in the surface mesh. The indices are the indices of the vertices in the surface mesh.
         * size: < spatial dim x number of vertices >
         */
        Matrix< DDRMat > mVertexCoordinates = Matrix< DDRMat >( 0, 0 );

        /**
         * @brief Displacements of the surface mesh vertices <dimension> x <number of vertices>
         *
         */
        Matrix< DDRMat > mDisplacements;

        /**
         * @brief List of vertices that are part of a given facet. The indices are the indices of the vertices in the surface mesh
         * size: number of facets< spatial_dim >
         */
        Vector< Vector< moris_index > > mFacetToVertexConnectivity;

        Vector< Vector< moris_index > > mVertexToVertexConnectivity;    // Input: vertex index, Output: All vertices connected by an edge to this vertex brendan documentation

        Vector< Vector< moris_index > > mVertexToFacetConnectivity;

        Matrix< DDRMat > mdShapeDiameterdVertex;    // cached shape diameter sensitivities


        /**
         * @brief Stores the facet normals for each facet in the surface mesh. The indices are the indices of the facets in the surface mesh, not the global indices!
         * size: < spatial dim x number of facets >
         */
        Matrix< DDRMat > mFacetNormals = Matrix< DDRMat >( 0, 0 );

      protected:    // variables
#if MORIS_HAVE_ARBORX
        /**
         * @brief ArborX Bounding volume hierarchy. Used to preselect which facets to check for intersection with a ray.
         *
         */
        ArborX::BVH< MemorySpace > mBVH;
#endif

        real mIntersectionTolerance = 1e-8;    // tolerance for interfaces when raycasting with this surface mesh
    };

}    // namespace moris::mtk

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
#include "cl_MTK_Shape_Diameter_Data.hpp"

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

    using Agglomeration_Function             = real ( * )( const real );    // Pointer to agglomeration function that takes a nodal shape diameter and returns an agglomerated value
    using Agglomeration_Sensitivity_Function = real ( * )( const real );    // Pointer to agglomeration sensitivity function that takes a nodal shape diameter and returns the sensitivity of the function wrt to the shape diameter

    struct Agglomeration_Parameters
    {
        real mExp;      // Exponent controls sharpness of agglomeration
        real mRef;      // Reference value to cut off any values greater than this
        real mShift;    // Shifts the value by this much

        Agglomeration_Parameters(
                real aAgglomerationExponent  = 2.0,
                real aAgglomerationReference = 1.0,
                real aAgglomerationShift     = 0.0 )
                : mExp( aAgglomerationExponent )
                , mRef( aAgglomerationReference )
                , mShift( aAgglomerationShift )
        {
            MORIS_ASSERT( (uint)mExp % 2 == 0, "Agglomeration_Parameters - Exponent must be even to ensure violation value is positive" );
        }
    };

    enum class Shape_Diameter_Method
    {
        RAYCAST,              // Use a weighted average of all ray intersections to compute the shape diameter
        INSCRIBED_CIRCLE,     // Use the diameter of the inscribed circle as the shape diameter (only for 2D meshes)
        SHORTEST_DISTANCE,    // Only use the closest ray intersection to compute the shape diameter
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
         * @brief Computes all intersection distances of a ray with the surface mesh
         * This method utilizes ArborX to find ray facet intersections, and then computes the intersection locations for the ray.
         * The method will only cast a single ray, which may not be sufficient to determine the region if the ray hits an edge or another pathological case is detected
         * The region is determined by the number of intersections. Even number = outside, Odd number = inside.
         *
         * @param aPoint Ray origin point. Passed by value as it may be altered
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @param aWarning Flag if the ray hits an edge or another pathological case
         * @param aIgnoreWarnings If true, all rays that hit a warning will be ceased immediately. Otherwise, the ray is cast and the result is returned as normal.
         * @return Intersection_Vector Vector of pairs <uint, real> which correspond to the facet index and the distance to the intersection point.
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
                const Matrix< DDRMat >&   aOrigins,
                const Matrix< DDRMat >&   aDirections,
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
                const Matrix< DDRMat >&           aOrigins,
                const Vector< Matrix< DDRMat > >& aDirections,
                Vector< Vector< bool > >&         aWarnings,
                bool                              aIgnoreWarnings = true ) const;

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
         * Casts a cone of rays from each node, centered around the vertex normal, and gets the distances for each ray-facet intersection
         *
         * @param aConeAngle Angle of the cone in degrees
         * @param aNumPolarRays Number of rays in the polar direction (θ)
         * @param aNumAzimuthRays Number of rays in the azimuth direction (φ) 1 if 2D
         * @return Shape_Diameter_Data struct which holds the weights and the intersection distances for each ray at each vertex
         */
        Shape_Diameter_Data cast_shape_diameter_ray_cones( const real aConeAngle, const uint aNumPolarRays, uint aNumAzimuthRays = 1 ) const;

        /**
         * Computes the shape diameter for each surface mesh node
         * The shape diameter is computed by casting a cone of rays from each node and using a weighted sum of each ray's closest distance
         * Quantity is computed per facet
         *
         * @param aConeAngle Angle of the cone in degrees
         * @param aNumPolarRays Number of rays in the polar direction (θ)
         * @param aNumAzimuthRays Number of rays in the azimuth direction (φ) 1 if 2D
         * @param
         * @return Vector< real > Shape diameter for each facet
         */
        Vector< real > compute_raycast_shape_diameter(
                const Agglomeration_Parameters& aAgglomeration,
                real                            aConeAngle,
                uint                            aNumPolarRays,
                uint                            aNumAzimuthRays = 1 );

        /**
         * Computes the shape diameter for each surface mesh node
         * The shape diameter is computed by computing the diameter of the smallest inscribed circle at each point. This approximates the medial axis distance
         * Quantity is originally computed per vertex, and then the average of all vertices for each facet is taken to give a quantity defined on facets
         *
         * @param aNumCircles Takes the mean of this many smallest circles for each point
         * @param aAngle Only chords between points within the cone are considered. This is because the sensitivity of the diameter goes to infinity if the chord is orthogonal to the vertex normal
         *
         * @return Vector< real > Shape diameter for each facet
         */
        Vector< real > compute_inscribed_circle_shape_diameter(
                const Agglomeration_Parameters& aAgglomeration,
                real                            aAngle         = 120.0,
                real                            aRelativeChord = 0.25,
                uint                            aNumCircles    = 1 );

        /**
         * Computes the shape diameter for each surface mesh node
         * The shape diameter is computed by finding the shortest distance to each facet and taking the minimum
         * Quantity is originally computed per vertex, and then the average of all vertices for each facet is taken to give a quantity defined on facets
         *
         * @param aNumCircles Takes the mean of this many shortest distances for each point
         * @param aAngle Ensures the direction vector of the shortest distance is within a cone of this angle from the vertex normal
         *
         * @return Vector< real > Shape diameter for each facet
         */
        Vector< real > compute_shortest_distance_shape_diameter(
                const Agglomeration_Parameters& aAgglomeration,
                real                            aAngle        = 40.0,
                uint                            aNumDistances = 1 );

        /**
         * Takes the shape diameter values computed on the facets and integrates them over the surface to get a global shape diameter violation
         *
         * Requires that the facet shape diameters are computed and stored in mShapeDiameterValues
         * @return real Integrated shape diameter violation over the surface mesh
         * */
        real integrate_shape_diameter_over_surface();

        /**
         * Computes a global shape diameter violation by integrating the shape diameter values computed on the facets with a raycast method
         * Rays are cast from the facet centers and the shape diameter is computed by a weighted average of the ray distances, and then integrated over the surface to get a global value
         *
         * @param aAgglom Agglomeration parameters struct to control the agglomeration function applied to the shape diameter values
         * @param aConeAngle Angle of the cone in degrees
         * @param aNumPolarRays Number of rays in the polar direction (θ)
         * @param aNumAzimuthRays Number of rays in the azimuth direction (φ) 1 if 2D
         * @return real Integrated shape diameter violation over the surface mesh
         */
        real compute_global_shape_diameter_raycast(
                const Agglomeration_Parameters& aAgglom,
                real                            aConeAngle,
                uint                            aNumPolarRays,
                uint                            aNumAzimuthRays = 1 );

        /**
         * Computes a global shape diameter violation by integrating the shape diameter values computed on the facets with an inscribed circle method
         * The smallest inscribed circle is computed for each vertex, and the shape diameter is computed by a weighted average of these distances, and then integrated over the surface to get a global value
         *
         * @param aAgglom Agglomeration parameters struct to control the agglomeration function applied to the shape diameter values
         * @param aConeAngle Angle of the cone in degrees
         * @param aRelativeChord Relative chord length to be used in inscribed circle computation
         * @return real Integrated shape diameter violation over the surface mesh
         */
        real compute_global_shape_diameter_inscribed_circle(
                const Agglomeration_Parameters& aAgglom,
                real                            aAngle         = 120.0,
                real                            aRelativeChord = 0.25,
                uint                            aNumCircles    = 1 );

        /**
         * Computes a global shape diameter violation by integrating the shape diameter values computed on the facets with a shortest distance method
         * The closet point to another part of the surface mesh within a cone is found for each vertex, and the shape diameter is computed by a weighted average of these distances, and then integrated over the surface to get a global value
         *
         * @param aAgglom Agglomeration parameters struct to control the agglomeration function applied to the shape diameter values
         * @param aConeAngle Angle of the cone in degrees
         * @param aNumDistances Number of closest distances to take the average of for each facet
         * @return real Integrated shape diameter violation over the surface mesh
         */
        real compute_global_shape_diameter_shortest_distance(
                const Agglomeration_Parameters& aAgglom,
                real                            aAngle        = 40.0,
                uint                            aNumDistances = 1 );

        /**
         * Computes the centroids of a single facet in the surface mesh
         * Size: <spatial dim> x <1>
         *
         * @return Matrix< DDRMat > facet centroid (average of the vertex coordinates of the facet)
         */
        Matrix< DDRMat > compute_facet_centroid( const uint aFacetIndex ) const;

        /**
         * Computes the centroids of all the facets in the surface mesh
         * Size: <spatial dim> x <number of facets>
         *
         * @return Matrix< DDRMat > facet centroids (average of the vertex coordinates of each facet)
         */
        Matrix< DDRMat > compute_facet_centroids() const;

        /**
         * Computes the measure (area in 3D, length in 2D) of a facet in the surface mesh
         * Virtual as some child implementations may prefer to compute and store this value. This implementation computes it on the fly.
         *
         * @param aFacetIndex local index of the facet to compute measure of.
         * @return Vector< real > facet measures. Size: <number of facets> x <1>
         */
        real compute_facet_measure( const uint aFacetIndex ) const;

        /**
         * Computes the measure (area in 3D, length in 2D) of all facets in the surface mesh
         * Virtual as some child implementations may prefer to compute and store this value. This implementation computes it on the fly.
         *
         * @return Vector< real > facet measures. Size: <number of facets> x <1>
         */
        virtual Vector< real > compute_facet_measure() const;

        /**
         * Computes the normal vector for each vertex as the average of the normals of its facets
         * Virtual as some child implementations may prefer to compute and store this value. This implementation computes it on the fly.
         *
         * @return Matrix< DDRMat > vertex normals stored column-wise. Size: <spatial dim> x <number of vertices>
         */
        Matrix< DDRMat > compute_vertex_normals();

        //-------------------------------------------------------------------------------
        // Quantity of interest sensitivities
        // -------------------------------------------------------------------------------

        /**
         * @brief Gets the derivative of the facet measure wrt a vertex's coordinates
         *
         * @param aFacetIndex local index of the facet to get sensitivities of
         * @param aVertexIndex local index of the vertex to get sensitivities of
         * @param aRequireIsMember If true, the method will check if the vertex is part of the facet and throw an error if not. If false, the method will return a zero matrix if the vertex is not part of the facet
         * @return Matrix< DDRMat > dMeasure/dVertex. Size: <1> x <spatial dim>
         */
        virtual Matrix< DDRMat > compute_dfacet_measure_dvertex( const uint aFacetIndex, const uint aVertexIndex, bool aRequireIsMember = false ) const;

        /**
         * Computes the derivative of the facet centroid wrt the coordinates of a facet vertex
         *
         * @param aFacetIndex local index of the facet to get sensitivities of
         * @param aVertexIndex local index of the vertex to get sensitivities with respect to
         * @param aRequireIsMember If true, the method will check if the vertex is part of the facet and throw an error if not. If false, the method will return a zero matrix if the vertex is not part of the facet
         *
         * @return Matrix< DDRMat > dCentroid/dVertex jacobian. Size: <spatial dim> x <spatial dim>
         * Rows correspond to the components of the centroid vector, and columns correspond to the components of the vertex coordinates
         */
        Matrix< DDRMat > compute_dfacet_centroid_dvertex( const uint aFacetIndex, const uint aVertexIndex, bool aRequireIsMember = false ) const;

        /**
         * Computes the derivative of the facet normal wrt the coordinates of a facet vertex
         *
         * @param aFacetIndex local index of the facet to get sensitivities of
         * @param aVertexIndex local index of the vertex to get sensitivities with respect to
         * @param aRequireIsMember If true, the method will check if the vertex is part of the facet and throw an error if not. If false, the method will return a zero matrix if the vertex is not part of the facet
         *
         * @return Matrix< DDRMat > dNormal/dVertex jacobian. Size: <spatial dim> x <spatial dim>
         * Rows correspond to the components of the normal vector, and columns correspond to the components of the vertex coordinates
         */
        Matrix< DDRMat > compute_dfacet_normal_dvertex( const uint aFacetIndex, const uint aVertexIndex, bool aRequireIsMember = false ) const;

        /**
         * Computes the derivative of the vertex normal wrt the coordinates of a vertex
         *
         * @param aVertexNormalIndex local index of the vertex to get the normal of
         * @param aVertexIndex local index of the vertex to get sensitivities with respect to. Each vertex normal depends on all vertices that share a facet with the vertex, so this function will return non-zero sensitivities for all of those vertices.
         * @param aRequireIsMember If true, the method will check if the vertex is part of any facet and throw an error if not. If false, the method will return a zero matrix if the vertex is not part of any facet
         *
         * @return Matrix< DDRMat > dNormal/dVertex jacobian. Size: <spatial dim> x <spatial dim>
         * Rows correspond to the components of the normal vector, and columns correspond to the components of the vertex coordinates
         */
        Matrix< DDRMat > compute_dvertex_normal_dvertex( const uint aVertexNormalIndex, uint aVertexIndex, bool aRequireIsMember = false ) const;

        /**
         * Computes the sensitivity of a ray-facet intersection distance with respect to the ray origin point
         *
         * @param aOrigin Origin point of the ray
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @param aFacetIndex Local index of the facet that the ray intersects. Must be a valid intersection
         * @return Matrix< DDRMat > dDistance/dOrigin gradient. Size: <1> x <spatial dim>
         */
        Matrix< DDRMat > compute_draycast_dorigin(
                const Matrix< DDRMat >& aOrigin,
                const Matrix< DDRMat >& aDirection,
                uint                    aFacetIndex ) const;

        /**
         * Computes the sensitivity of a ray-facet intersection distance with respect to the ray direction vector
         *
         * @param aOrigin Origin point of the ray
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @param aFacetIndex Local index of the facet that the ray intersects. Must be a valid intersection
         * @return Matrix< DDRMat > dDistance/dDirection gradient. Size: <1> x <spatial dim>
         */
        Matrix< DDRMat > compute_draycast_ddirection(
                const Matrix< DDRMat >& aOrigin,
                const Matrix< DDRMat >& aDirection,
                uint                    aFacetIndex ) const;

        /**
         * Computes the sensitivity of a ray-facet intersection distance with respect to ALL vertices of the facet
         *
         * @param aOrigin Origin point of the ray
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @param aFacetIndex Local index of the facet that the ray intersects. Must be a valid intersection
         * @return Matrix< DDRMat > dDistance/dOrigin gradient. Size: <spatial dim> x <spatial dim>
         * Columns correspond to the components of the vertex coordinates, rows correspond to the vertices of the facet
         */
        Matrix< DDRMat > compute_draycast_dvertices(
                const Matrix< DDRMat >& aOrigin,
                const Matrix< DDRMat >& aDirection,
                uint                    aFacetIndex ) const;

        /**
         * Gets the sum of the nodal shape diameter sensitivities wrt to all vertex coordinates
         *
         * @return Matrix< DDRMat > dDiameter/dVertices. Size: < number of vertices > x <spatial dim>
         */
        const Matrix< DDRMat >& get_facet_shape_diameter_sensitivities() const;

        /**
         * Gets the shape diameter of a single facet in the surface mesh
         *
         * @param aFacetIndex local index of the facet
         * @return real shape diameter of the facet
         */
        const real get_facet_shape_diameter( uint aFacetIndex ) const;

        /**
         * Gets the derivative of the shape diameter wrt a vertex's coordinates
         *
         * @param aVertexIndex local index of the vertex to get sensitivities of
         * @return Matrix< DDRMat > dDiameter/dVertex. Size: <1> x <spatial dim>
         */
        Matrix< DDRMat > compute_ddiameter_dvertex( uint aVertexIndex ) const;

        static real tanh_clip( real aShapeDiameter, const Agglomeration_Parameters& aAgglom );

        static real dtanh_clip( real aShapeDiameter, const Agglomeration_Parameters& aAgglom );

        static real max_clip( real aShapeDiameter, const Agglomeration_Parameters& aAgglom );

        static real dmax_clip( real aShapeDiameter, const Agglomeration_Parameters& aAgglom );

        static real sech( real aValue );

        static real normal( real aValue, real aStdDev );    // assumes 0 mean

        static real dnormal_dsigma( real aValue, real aStdDev );    // assumes 0 mean

        /**
         * @brief Computes the absolute value of the angle between two vectors in radians
         * @param aVectorA First vector
         * @param aVectorB Second vector
         * @return Angle between the two vectors in radians
         */
        static real compute_dot_absolute(
                const Matrix< DDRMat >& aVectorA,
                const Matrix< DDRMat >& aVectorB );

        /**
         * @brief Computes the derivative of the dot product between two vectors with respect to the first vector
         */
        static Matrix< DDRMat > compute_ddot_absolute(
                const Matrix< DDRMat >& aVectorA,
                const Matrix< DDRMat >& aVectorB );

        /**
         * @brief Takes a quantity defined on vertices and averages it over the facets to give a quantity defined on facets
         *
         * @tparam T Type of the quantity to be averaged. Must have += and / operators defined.
         * @param aVertexQuantity Vector of the quantity defined on vertices. Size: <number of vertices>
         * @return Vector< real > quantity averaged on facets. Size: <number of facets>
         */
        template< typename T >
        Vector< T > average_vertex_quantity_over_facets( const Vector< T >& aVertexQuantity ) const
        {
            MORIS_ERROR( aVertexQuantity.size() == this->get_number_of_vertices(),
                    "Surface_Mesh::average_vertex_quantity_over_facets - The input vertex quantity size %ld does not match the number of vertices in the mesh %d.",
                    aVertexQuantity.size(),
                    this->get_number_of_vertices() );

            Vector< T > tFacetQuantity( this->get_number_of_facets() );

            for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
            {
                const Vector< moris_index >& tFacetVertices = this->get_facets_vertex_indices( iF );
                T                            tSum           = 0.0;
                for ( uint iV : tFacetVertices )
                {
                    tSum += aVertexQuantity( iV );
                }
                tFacetQuantity( iF ) = tSum / (real)tFacetVertices.size();
            }

            return tFacetQuantity;
        }

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
         * @brief computes the normal vector of a single facet in the surface mesh. Stores the result in mFacetNormals.
         * Using the return value of this function is the same as calling get_facet_normal( aFacetIndex ) after this function
         *
         * @param aFacetIndex local facet index to compute normal for
         */
        Matrix< DDRMat > compute_facet_normal( const uint aFacetIndex );

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

        // static real dstdev_ddata( real aMean );

        // static real dmean_ddata( uint aNumData );

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
                ExecutionSpace const &            aExecutionSpace,
                const Matrix< DDRMat >&           aOrigins,
                const Vector< Matrix< DDRMat > >& aDirections );

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
                ExecutionSpace const &  aExecutionSpace,
                const Matrix< DDRMat >& aOrigins,
                const Matrix< DDRMat >& aDirections );

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
                const Matrix< DDRMat >&           aOrigins,
                const Vector< Matrix< DDRMat > >& aDirections ) const;

        /**
         * @brief Uses ArborX bounding volume hierarchy to determine which facets may be intersected by the ray
         * The ray only travels in the positive direction.
         *
         * @param aOrigins Origin points for the rays. Each column is a new origin, size <dimension> x <number of origins>
         * @param aDirections Directions for each ray. This version casts the same directions for every origin.size <dimension> x <number of directions>
         * @return Vector< Vector< Vector< uint > > > Innermost vector is the indices of the facets that the ray hit. Middle vector is for each direction for each origin point. Outer vector is for origin points
         */
        Vector< Vector< Vector< uint > > > batch_preselect_with_arborx(
                const Matrix< DDRMat >& aOrigins,
                const Matrix< DDRMat >& aDirections ) const;

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

        Vector< Vector< moris_index > > mVertexToVertexConnectivity;    // Input: vertex index, Output: All vertices connected by an edge to this vertex

        Vector< Vector< moris_index > > mVertexToFacetConnectivity;

        real             mIntegratedShapeDiameter = MORIS_REAL_MAX;    // cached global shape diameter value
        Vector< real >   mShapeDiameters;                              // cached nodal shape diameter values
        Matrix< DDRMat > mdShapeDiameterdVertex;                       // cached shape diameter sensitivities

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

        mutable real mIntersectionTolerance = 1e-8;    // tolerance for interfaces when raycasting with this surface mesh
    };

}    // namespace moris::mtk

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
#include "fn_MTK_QuadraturePointMapper_Ray_ArborX_Details.hpp"

#include <ArborX.hpp>
#include <ArborX_Box.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_View.hpp>
#include <decl/Kokkos_Declare_SERIAL.hpp>

namespace moris::mtk
{
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using MemorySpace    = ExecutionSpace::memory_space;

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
        Surface_Mesh( Matrix< DDRMat >          aVertexCoordinates,
                Vector< Vector< moris_index > > aFacetConnectivity,
                real                            aIntersectionTolerance );

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
         * Adds a displacement to the current displacement of the vertex
         *
         * @param aVertexIndex vertex index to add displacement to
         * @param aDisplacement <dimension > x < 1 > matrix containing displacment to add
         */
        void append_vertex_displacement( const uint aVertexIndex, const Matrix< DDRMat >& aDisplacement );

        /**
         * @brief removes any displacement, rotation, and scaling applied to the surface mesh
         */
        void reset_coordinates();

        // -------------------------------------------------------------------------------
        // Accessor methods
        // -------------------------------------------------------------------------------

        [[nodiscard]] virtual Matrix< DDRMat > get_all_vertex_coordinates() const;

        /**
         * @brief Gets the coordinates of a single vertex from the local index aVertexIndex
         */
        [[nodiscard]] virtual Matrix< DDRMat > get_vertex_coordinates( const uint aVertexIndex ) const;

        /**
         * @brief Gets the indices to the vertices that form the facet with the local index aFacetIndex
         *
         * @param aFacetIndex local index of the facet
         * @return Vector< moris_index > local vertex indices that form the facet
         */
        [[nodiscard]] Vector< moris_index > get_facets_vertex_indices( const uint aFacetIndex ) const;

        /**
         * @brief Gets the coordinates of all vertices that form the facet with the local index aFacetIndex
         * Size: < spatial dim x number of vertices in the facet >
         */
        [[nodiscard]] Matrix< DDRMat > get_all_vertex_coordinates_of_facet( const uint aFacetIndex ) const;

        /**
         * @brief Returns the facet normals for each facet in the surface mesh.
         * @return A (d x n) matrix where d is the dimension of the mesh and n is the number of facets in the surface mesh (holding the normal components).
         */
        [[nodiscard]] Matrix< DDRMat > get_all_facet_normals() const;

        /**
         * @brief Gets the normal vector of the facet with the local index aFacetIndex
         *
         * @param aFacetIndex local index of the facet
         * @return Matrix< DDRMat > normal vector of the facet
         */
        [[nodiscard]] Matrix< DDRMat > get_facet_normal( const uint aFacetIndex ) const;

        [[nodiscard]] virtual uint get_spatial_dimension() const;

        /**
         * @brief Returns the number of facets in the mesh
         */
        [[nodiscard]] uint get_number_of_facets() const;

        [[nodiscard]] uint get_number_of_vertices() const;


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
        raycast_point( const Matrix< DDRMat >& aPoint ) const;

        /**
         * @brief Gets the intersection distances of all the facets intersected by the infinite ray
         *
         * @param aPoint origin point of the ray
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @param aIntersectionFacetIndices Return value. Indices of the facets that the ray hit
         * @return Each index of the output is a distance along the ray that a facet intersects with.
         * The global coordinates of the intersection can be found by computing aPoint + <output_entry> * aDirection for each entry of the vector
         *
         */
        Vector< real > compute_ray_facet_intersections(
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection,
                Vector< uint >&         aIntersectionFacetIndices ) const;

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
                const Matrix< DDRMat >& aDirection ) const;

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
                const Matrix< DDRMat >& aDirection ) const;

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
                const Matrix< DDRMat >& aDirection ) const;

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


        //-------------------------------------------------------------------------------
        // Initialization methods
        // -------------------------------------------------------------------------------

      protected:    // methods
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
         * @brief Constructs the ArborX bounding volume hierarchy for the surface mesh and stores it in mBVH
         *
         */
        void construct_bvh();

      private:
        /**
         * @brief Determines if a point is inside or outside the surface mesh via raycasting
         * This method utilizes ArborX to find ray facet intersections, and then computes the intersection locations for the ray.
         * The method will only cast a single ray, which may not be sufficient to determine the region if the ray hits an edge or another pathological case is detected
         * The region is determined by the number of intersections. Even number = outside, Odd number = inside.
         *
         * @param aPoint Ray origin point. Passed by value as it may be altered
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @param aIntersectionFacetIndices Return value. Indices of the facets that the ray hit
         */
        Vector< real >
        cast_single_ray(
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection,
                Vector< uint >&         aIntersectionFacetIndices ) const;

        /**
         * @brief Determines if a point is inside or outside the surface mesh via raycasting
         * This method utilizes ArborX to find ray facet intersections, and then computes the intersection locations for the ray.
         * The method will only cast a single ray, which may not be sufficient to determine the region if the ray hits an edge or another pathological case is detected
         * The region is determined by the number of intersections. Even number = outside, Odd number = inside.
         *
         * @param aPoint Ray origin point. Passed by value as it may be altered
         * @param aDirection Direction that the ray casts in. Does not have to be a unit vector
         * @param aIntersectionFacetIndices Return value. Indices of the facets that the ray hit
         */
        Vector< real >
        cast_single_ray(
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aDirection ) const;

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

      protected:    // variables
        /**
         * @brief List of cell indices that the vertex with the given index is part of. The indices are the indices of
         * the cell in the surface mesh
         * size: number of facets< spatial_dim >
         */
        Vector< Vector< moris_index > > mFacetConnectivity;

        /**
         * @brief ArborX Bounding volume hierarchy. Used to preselect which facets to check for intersection with a ray.
         *
         */
        ArborX::BVH< MemorySpace > mBVH;

        /**
         * @brief Stores the facet normals for each facet in the surface mesh. The indices are the indices of the facets in the surface mesh, not the global indices!
         * size: < spatial dim x number of facets >
         */
        Matrix< DDRMat > mFacetNormals = Matrix< DDRMat >( 0, 0 );

        real mIntersectionTolerance = 1e-8;    // tolerance for interfaces when raycasting with this surface mesh
    };

}    // namespace moris::mtk
/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Surface_Mesh_QI.hpp
 *
 */
#pragma once

#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Shape_Diameter_Data.hpp"
#include "cl_MTK_Agglomeration_Parameters.hpp"

namespace moris::sol
{
    class Dist_Vector;
}

namespace moris::mtk
{
    enum class Shape_Diameter_Method
    {
        RAYCAST,              // Use a weighted average of all ray intersections to compute the shape diameter
        INSCRIBED_CIRCLE,     // Use the diameter of the inscribed circle as the shape diameter (only for 2D meshes)
        SHORTEST_DISTANCE,    // Only use the closest ray intersection to compute the shape diameter
    };

    // User callbacks receive the complete, read-only surface mesh API.
    using Surface_Mesh_QI_Function             = real ( * )( const Surface_Mesh& );
    using Surface_Mesh_QI_Sensitivity_Function = Matrix< DDRMat > ( * )( const Surface_Mesh&, uint );

    /**
     * Computes surface mesh quantities of interest and owns their intermediate values.
     * The referenced mesh must outlive this object. Recompute the QI after changing the mesh
     * before requesting its sensitivities. Separate QIs can share a mesh without sharing caches.
     */
    class Surface_Mesh_QI
    {
      public:
        explicit Surface_Mesh_QI( const Surface_Mesh& aSurfaceMesh );

        /**
         * Selects and evaluates a QI, loading user callbacks from the input library if needed.
         * USER_DEFINED requires function_name and sensitivity_function_name. The value callback
         * returns a real; its sensitivity callback returns a 1 x spatial dimension row matrix
         * for the supplied local surface vertex index. Callbacks may cache data in the input file.
         */
        real compute_QI( QI_Type aType, const Parameter_List& aParameters, const std::shared_ptr< Library_IO >& aLibrary );

        /**
         * Evaluates callbacks supplied directly, using the same workflow as loaded callbacks.
         */
        real compute_QI( Surface_Mesh_QI_Function    aFunction,
                Surface_Mesh_QI_Sensitivity_Function aSensitivityFunction );

        /**
         * Returns the sensitivity of the last compute_QI call as a 1 x spatial dimension row matrix.
         */
        Matrix< DDRMat > compute_dQI_dvertex( uint aVertexIndex ) const;

        /**
         * Adds QI sensitivities to the requested distributed vector column.
         * PDV IDs are indexed by spatial dimension, then local surface vertex. -1 skips a PDV.
         */
        void compute_QI_sensitivities( const Vector< Vector< moris_index > >& aVertexPDVIDs,
                sol::Dist_Vector*                                             aSensitivities,
                uint                                                          aRequestIndex ) const;

        /**
         * Computes the volume enclosed by the surface mesh
         */
        real compute_volume() const;

        /**
         * Approximates the area between a selected 2D source loop and a fixed target boundary.
         * Rays follow the negative length-weighted vertex normals; source facet orientation
         * must make these directions point toward the target. Each source facet and its two
         * nearest positive ray intersections define a quadrilateral whose absolute shoelace
         * area is accumulated. Source vertices share their ray hits across incident facets.
         *
         * @param aSourceFacets Local facet indices of the moving source loop.
         * @param aVertexSensitivities Output dArea/dVertex, size <number of vertices> x <2>.
         * Only source-loop coordinates vary; target-boundary rows remain zero. Includes direct
         * vertex motion, ray-origin motion, ray-direction motion, and changing normal weights.
         * The derivative applies while the intersected facets and quadrilateral orientations stay fixed.
         */
        real compute_raycast_wedge_area( const Vector< moris_index >& aSourceFacets,
                Matrix< DDRMat >&                                     aVertexSensitivities ) const;

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

        static real tanh_clip( real aShapeDiameter, const Agglomeration_Parameters& aAgglom, const Matrix< DDRMat >& aCoords, uint aIndex = 0 );

        static real dtanh_clip( real aShapeDiameter, const Agglomeration_Parameters& aAgglom, const Matrix< DDRMat >& aCoords, uint aIndex = 0 );

        static real max_clip( real aShapeDiameter, const Agglomeration_Parameters& aAgglom, const Matrix< DDRMat >& aCoords, uint aIndex = 0 );

        static real dmax_clip( real aShapeDiameter, const Agglomeration_Parameters& aAgglom, const Matrix< DDRMat >& aCoords, uint aIndex = 0 );

        static real sech( real aValue );

        static real dnormal_dsigma( real aValue, real aStdDev );    // assumes 0 mean

        void write_shape_diameters_to_file( const std::string& aFilePath ) const;

      private:
        const Surface_Mesh&                  mSurfaceMesh;
        QI_Type                              mQIType              = QI_Type::END_ENUM;
        Surface_Mesh_QI_Sensitivity_Function mSensitivityFunction = nullptr;

        real             mIntegratedShapeDiameter = MORIS_REAL_MAX;
        Vector< real >   mShapeDiameters;
        Matrix< DDRMat > mdShapeDiameterdVertex;
    };
}    // namespace moris::mtk

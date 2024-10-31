/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Raycaster.cpp
 *
 */

// #pragma once

// #include <fstream>

// #include "cl_Matrix.hpp"
// #include "linalg_typedefs.hpp"
// #include "moris_typedefs.hpp"
// #include "SDF_Tools.hpp"

// #include "cl_MTK_Surface_Mesh.hpp"

// namespace moris::sdf
// {
//     /**
//      * Conducts a raycast to determine if aPoint is inside or outside of the object. Casts rays in each coordinate direction until the point is resolved.
//      * Also determines intersection locations along the coordinate axes of aPoint with each facet.
//      *
//      * @param aPoint coordinate point in which the ray will originate. Passed by value to avoid modifying the original point if rotations are necessary
//      */
//     Object_Region
//     get_region_from_raycast(
//             mtk::Surface_Mesh& aMesh,
//             Matrix< DDRMat >   aPoint );

//-------------------------------------------------------------------------------

//     /**
//      * Casts in the aAxis direction and returns the distance to every facet the ray hits
//      *
//      * @param aMesh water tight collection of facets to cast on to
//      * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
//      * @param aAxis direction in which the ray is cast
//      * @param aFacetIndices return value. The facet indices associated with the intersection locations
//      * @return intersection coordinate locations
//      */
//     Vector< real >
//     compute_distance_to_facets(
//             mtk::Surface_Mesh&      aMesh,
//             const Matrix< DDRMat >& aPoint,
//             const uint              aAxis,
//             Vector< uint >&         aIntersectedFacets );

//-------------------------------------------------------------------------------

//     /**
//      * Performs a raycast in the aAxis direction and determines whether the point is inside or outside.
//      * Can be used for 2D or 3D geometries
//      *
//      * @param aMesh water tight collection of facets to cast on to
//      * @param aFacetMinCoords minimum bounding coordinates of aMesh
//      * @param aFacetMaxCoords maximum bounding coordinates of aMesh
//      * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
//      * @param aAxis direction in which the ray is cast
//      * @return whether the point is inside aMesh or not. 0 = outside, 1 = inside, 2 = unsure
//      */
//     Object_Region
//     voxelize(
//             mtk::Surface_Mesh& aMesh,
//             Matrix< DDRMat >&  aPoint,
//             uint               aAxis );

//-------------------------------------------------------------------------------

//     /**
//      * @brief Takes a point in 3D space and determines which triangles
//      * in the positive and negative x direction the point could possibly intersect.
//      * These triangles are returned. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
//      *
//      * @param aMesh water tight collection of facets to cast on to
//      * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
//      * @param aAxis direction in which the ray is cast
//      * @param aCandidateFacets return variable, facet indices which aPoint lies in the bounding box of, and intersection coordinates need to be computed to determine intersection
//      * @return Detection of pathological preselection cases, like a ray cast directly onto a vertex or the cast point being on a vertex
//      */
//     Preselection_Result
//     preselect_triangles(
//             mtk::Surface_Mesh&      aMesh,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis,
//             Vector< uint >&         aCandidateFacets );


//     /**
//      * @brief Takes a point in 3D space and determines which triangles
//      * in the positive and negative x direction the point could possibly intersect.
//      * These triangles are returned. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
//      *
//      * @param aMesh water tight collection of facets to cast on to
//      * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
//      * @param aAxis direction in which the ray is cast
//      * @param aCandidateFacets return variable, facet indices which aPoint lies in the bounding box of, and intersection coordinates need to be computed to determine intersection
//      * @return Detection of pathological preselection cases, like a ray cast directly onto a vertex or the cast point being on a vertex
//      */
//     void
//     preselect_triangles_moller_trumbore(
//             mtk::Surface_Mesh&      aMesh,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis,
//             Vector< uint >&         aCandidateFacets );

//-------------------------------------------------------------------------------

//     /**
//      * Determines which lines are definitely pierced by the ray and which lines could be pierced by the ray. Returns true only if all preselection was successful
//      *
//      * @param aMesh water tight collection of facets to cast on to
//      * @param aFacetMinCoords minimum bounding coordinates of aMesh
//      * @param aFacetMaxCoords maximum bounding coordinates of aMesh
//      * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
//      * @param aAxis direction in which the ray is cast
//      * @param aIntersectedFacets return variable, facets which are for sure intersected by the ray
//      * @param aCandidateFacets return variable, facets which aPoint lies in the bounding box of, and intersection coordinates need to be computed to determine intersection
//      * @return Detection of pathological preselection cases, like a ray cast directly onto a vertex or the cast point being on a vertex
//      */
//     Preselection_Result
//     preselect_lines(
//             mtk::Surface_Mesh&      aMesh,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis,
//             Vector< uint >&         aIntersectedFacets,
//             Vector< uint >&         aCandidateFacets );

//-------------------------------------------------------------------------------

//     /**
//      * checks if the cast point will cast a ray that will actually intersect aCandidateFacets. If so, the facet is added to the return
//      *
//      * @param aCandidateFacets indices of facets which to determine intersection
//      * @param aMesh water tight collection of facets to cast on to
//      * @param aPoint spatial location of the origin of the ray
//      * @param aAxis coordinate axis in which to cast from
//      * @return pointers to which facets are intersected
//      */
//     Vector< uint >
//     intersect_triangles(
//             Vector< uint >&         aCandidateFacets,
//             mtk::Surface_Mesh&                 aMesh,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis );

//-------------------------------------------------------------------------------

//     /**
//      * checks if the cast point will cast a ray that will actually intersect aCandidateFacets. If so, the facet is added to the return
//      *
//      * @param aCandidateFacets indices of facets which to determine intersection
//      * @param aMesh water tight collection of facets to cast on to
//      * @param aPoint spatial location of the origin of the ray
//      * @param aAxis coordinate axis in which to cast from
//      * @param aIntersectedFacets return variable, facets which are for sure intersected by the ray
//      */
//     Vector< real >
//     intersect_lines_moller_trumbore(
//             Vector< uint >&         aCandidateFacets,
//             mtk::Surface_Mesh&      aMesh,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis,
//             Vector< uint >&         aIntersectedFacets );


//-------------------------------------------------------------------------------

//     /**
//      * checks if the cast point will cast a ray that will actually intersect aCandidateFacets. If so, the facet is added to the return
//      *
//      * @param aCandidateFacets indices of facets which to determine intersection
//      * @param aMesh water tight collection of facets to cast on to
//      * @param aPoint spatial location of the origin of the ray
//      * @param aAxis coordinate axis in which to cast from
//      * @param aIntersectedFacets return variable, facets which are for sure intersected by the ray
//      */
//     Vector< real >
//     intersect_triangles_moller_trumbore(
//             Vector< uint >&         aCandidateFacets,
//             mtk::Surface_Mesh&      aMesh,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis,
//             Vector< uint >&         aIntersectedFacets );

//-------------------------------------------------------------------------------

//     /**
//      * Takes all of potential facets (in aIntersectedFacets) and computes the coordinate axis intersection location
//      * with the ray originating from aPoint. Does not handle cases where the ray intersects with the vertex
//      *
//      * @param aIntersectedFacets facets which are for sure intersected by the ray, but the intersection coordinate is not known
//      * @param aPoint spatial location of the origin of the ray
//      * @param aAxis coordinate axis to shoot ray in. 0 = x, 1 = y, 2 = z
//      * @return matrix of intersection coordinate locations
//      */
//     Vector< real >
//     intersect_ray_with_facets(
//             mtk::Surface_Mesh&      aMesh,
//             Vector< uint >&         aIntersectedFacets,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis,
//             Preselection_Result     aPreselection,
//             bool                    aIgnoreErrors = false );

//-------------------------------------------------------------------------------

//     /**
//      * determines if aPoint is within the closed collection of triangles
//      * NOTE: Only to be used with 3D raycasts
//      *
//      * @param aIntersectionCoords locations where the ray intersects the facets
//      * @param aPoint spatial location of the origin of the ray
//      * @param aAxis coordinate axis in which the ray was cast from
//      * @return whether the point is inside the object or not. 0 = outside, 1 = inside, 2 = unsure
//      */
//     Object_Region
//     check_if_node_is_inside_triangles(
//             Vector< real >&         aIntersectionCoords,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis,
//             const real&             aIntersectionTolerance = 1e-8 );

//-------------------------------------------------------------------------------

//     /**
//      * determines if aPoint is within the closed collection of lines.
//      * NOTE: Only to be used with 2D raycasts
//      * NOTE: aIntersectionCoords and aCandidateFacets contain mutually exclusive lists of facets.
//      *
//      * @param aIntersectionCoords locations where the ray intersects the facets
//      * @param aCandidateFacets facets which are for sure intersected by the ray
//      * @param aPoint spatial location of the origin of the ray
//      * @param aAxis coordinate axis in which the ray was cast from
//      * @return whether the point is inside the object or not. 0 = outside, 1 = inside, 2 = unsure
//      */
//     Object_Region
//     check_if_node_is_inside_lines(
//             mtk::Surface_Mesh&      aMesh,
//             const Vector< real >&   aIntersectionCoords,
//             const Vector< uint >&   aCandidateFacets,
//             const Matrix< DDRMat >& aPoint,
//             uint                    aAxis );

//-------------------------------------------------------------------------------

//     /**
//      * Generates a random rotation angle (and axis in 3D) and builds a rotation matrix.
//      * Rotates both aMesh and aPoint to raycast in another random direction in the event the node is unsure.
//      *
//      * @param aMesh the object to rotate
//      * @param aPoint the cast point to rotate
//      */
//     void
//     random_rotation(
//             mtk::Surface_Mesh& aMesh,
//             Matrix< DDRMat >&  aPoint );

//-------------------------------------------------------------------------------

// }    // namespace moris::sdf

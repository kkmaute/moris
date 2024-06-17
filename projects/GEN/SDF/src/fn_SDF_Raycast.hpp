/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Raycaster.cpp
 *
 */

#pragma once

#include <fstream>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "SDF_Tools.hpp"

#include "cl_SDF_Object.hpp"

namespace moris::sdf
{
    /**
     * Conducts a raycast to determine if aPoint is inside or outside of the object. Casts rays in each coordinate direction until the point is resolved.
     * Also determines intersection locations along the coordinate axes of aPoint with each facet.
     *
     * @param aPoint coordinate point in which the ray will originate
     */
    Object_Region
    raycast_point(
            Object&          aObject,
            Matrix< DDRMat > aPoint );

    //-------------------------------------------------------------------------------

    /**
     * Casts in the aAxis direction and returns the distance to every facet the ray hits
     *
     * @param aObject water tight collection of facets to cast on to
     * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
     * @param aAxis direction in which the ray is cast
     * @param aFacetIndices return value. The facet indices associated with the intersection locations
     * @return intersection coordinate locations
     */
    Vector< real >
    compute_distance_to_facets(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            const uint              aAxis,
            Vector< Facet* >&       aIntersectedFacets );

    //-------------------------------------------------------------------------------

    /**
     * Performs a raycast in the aAxis direction and determines whether the point is inside or outside.
     * Can be used for 2D or 3D geometries
     *
     * @param aObject water tight collection of facets to cast on to
     * @param aFacetMinCoords minimum bounding coordinates of aObject
     * @param aFacetMaxCoords maximum bounding coordinates of aObject
     * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
     * @param aAxis direction in which the ray is cast
     * @return whether the point is inside aObject or not. 0 = outside, 1 = inside, 2 = unsure
     */
    Object_Region
    voxelize(
            Object&           aObject,
            Matrix< DDRMat >& aPoint,
            uint              aAxis );

    //-------------------------------------------------------------------------------

    /**
     * @brief Takes a point in 3D space and determines which triangles
     * in the positive and negative x direction the point could possibly intersect.
     * These triangles are returned. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
     *
     * @param aObject water tight collection of facets to cast on to
     * @param aFacetMinCoords minimum bounding coordinates of aObject
     * @param aFacetMaxCoords maximum bounding coordinates of aObject
     * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
     * @param aAxis direction in which the ray is cast
     * @return facet indices which could be intersected by the ray
     */
    Vector< uint >
    preselect_triangles(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis );

    //-------------------------------------------------------------------------------

    /**
     * Determines which lines are definitely pierced by the ray and which lines could be pierced by the ray. Returns true only if all preselection was successful
     *
     * @param aObject water tight collection of facets to cast on to
     * @param aFacetMinCoords minimum bounding coordinates of aObject
     * @param aFacetMaxCoords maximum bounding coordinates of aObject
     * @param aPoint Point in space that lies within the bounding plane of the candidate triangles
     * @param aAxis direction in which the ray is cast
     * @param aIntersectedFacets return variable, facets which are for sure intersected by the ray
     * @param aCandidateFacets return variable, facets which aPoint lies in the bounding box of, and intersection coordinates need to be computed to determine intersection
     * @return true if the preselection did not find any facets whose vertex would be hit by cast, false otherwise
     */
    Preselection_Result
    preselect_lines(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            Vector< uint >&         aIntersectedFacets,
            Vector< Facet* >&       aCandidateFacets );

    //-------------------------------------------------------------------------------


    /**
     * checks if the cast point will cast a ray that will actually intersect aCandidateFacets. If so, the facet is added to the return
     *
     * @param aCandidateFacets indices of facets which to determine intersection
     * @param aObject water tight collection of facets to cast on to
     * @param aPoint spatial location of the origin of the ray
     * @param aAxis coordinate axis in which to cast from
     * @return pointers to which facets are intersected
     */
    Vector< Facet* >
    intersect_triangles(
            Vector< uint >&         aCandidateFacets,
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis );

    //-------------------------------------------------------------------------------

    /**
     * Takes all of potential facets (in aIntersectedFacets) and computes the coordinate axis intersection location
     * with the ray originating from aPoint. Does not handle cases where the ray intersects with the vertex
     *
     * @param aIntersectedFacets facets which are for sure intersected by the ray, but the intersection coordinate is not known
     * @param aPoint spatial location of the origin of the ray
     * @param aAxis coordinate axis to shoot ray in. 0 = x, 1 = y, 2 = z
     * @return matrix of intersection coordinate locations
     */
    Vector< real >
    intersect_ray_with_facets(
            Vector< Facet* >&       aIntersectedFacets,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            Preselection_Result     aRayOnVertex = SUCCESS );

    //-------------------------------------------------------------------------------

    /**
     * determines if aPoint is within the closed collection of triangles
     * NOTE: Only to be used with 3D raycasts
     *
     * @param aIntersectionCoords locations where the ray intersects the facets
     * @param aPoint spatial location of the origin of the ray
     * @param aAxis coordinate axis in which the ray was cast from
     * @return whether the point is inside the object or not. 0 = outside, 1 = inside, 2 = unsure
     */
    Object_Region
    check_if_node_is_inside_triangles(
            Vector< real >&   aIntersectionCoords,
            Matrix< DDRMat >& aPoint,
            uint              aAxis );

    //-------------------------------------------------------------------------------

    /**
     * determines if aPoint is within the closed collection of lines.
     * NOTE: Only to be used with 2D raycasts
     * NOTE: aIntersectionCoords and aCandidateFacets contain mutually exclusive lists of facets.
     *
     * @param aIntersectionCoords locations where the ray intersects the facets
     * @param aCandidateFacets facets which are for sure intersected by the ray
     * @param aPoint spatial location of the origin of the ray
     * @param aAxis coordinate axis in which the ray was cast from
     * @return whether the point is inside the object or not. 0 = outside, 1 = inside, 2 = unsure
     */
    Object_Region
    check_if_node_is_inside_lines(
            Object&                 aObject,
            const Vector< real >&   aIntersectionCoords,
            const Vector< uint >&   aCandidateFacets,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis );

    //-------------------------------------------------------------------------------

    /**
     * Generates a random rotation angle (and axis in 3D) and builds a rotation matrix.
     * Rotates both aObject and aPoint to raycast in another random direction in the event the node is unsure.
     *
     * @param aObject the object to rotate
     * @param aPoint the cast point to rotate
     */
    void
    random_rotation(
            Object&           aObject,
            Matrix< DDRMat >& aPoint );

    //-------------------------------------------------------------------------------

}    // namespace moris::sdf
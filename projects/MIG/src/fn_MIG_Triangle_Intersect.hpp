/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_MIG_Triangle_Intersect.hpp
 *
 */

#pragma once
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "moris_typedefs.hpp"

namespace moris::mig
{
    /*
     * All intersection points of two triangles
     * @param[ in ]  aFirstTRICoords coordinates of first triangle
     * @param[ in ]  aSecondTRICoords coordinates of second triangle
     * @param[ out ] aIntersectedPoints intersection points size <2, num of intersections>
     */
    void
    Intersect(
            Matrix< DDRMat > const &aFirstTRICoords,
            Matrix< DDRMat > const &aSecondTRICoords,
            Matrix< DDRMat >       &aIntersectedPoints );

    /*
     * Checks if two triangles intersect
     * @param[ in ]  aFirstTRICoords coordinates of first triangle
     * @param[ in ]  aSecondTRICoords coordinates of second triangle
     * @return true if any part of the triangles intersect or one triangle is completely inside the other
     */
    bool triangles_intersect(
            Matrix< DDRMat > const &aFirstTRICoords,
            Matrix< DDRMat > const &aSecondTRICoords );

    /*
     * Computes locations of intersections along the edges of two triangles
     * @param[ in ]  aFirstTRICoords coordinates of first triangle
     * @param[ in ]  aSecondTRICoords coordinates of second triangle
     * @param[ out ] aIntersectedPoints intersection points size <2, num of intersections>
     */
    void
    edge_intersect(
            Matrix< DDRMat > const &aFirstTRICoords,
            Matrix< DDRMat > const &aSecondTRICoords,
            Matrix< DDRMat >       &aIntersectedPoints );

    /*
     * Finds vertices of one triangle within another one
     * @param[ in ]  aFirstTRICoords coordinates of triangle to check for insidedness
     * @param[ in ]  aSecondTRICoords coordinates of triangle to check against
     * @param[ out ] aIntersectedPoints Vertices of first triangle which are inside the second triangle
     * @note (point coordinates are stored column-wise, in counter clock
     *order) the corners of first which lie in the interior of second.
     */
    void find_vertices_inside_triangle(
            Matrix< DDRMat > const &aFirstTRICoords,
            Matrix< DDRMat > const &aSecondTRICoords,
            Matrix< DDRMat >       &aIntersectedPoints );

    /*
     * sort points and remove duplicates
     * orders polygon corners in counter clock wise and removes duplicates
     * @param[ in ] aIntersectedPoints polygon points, not ordered
     * @param[ out ] aIntersectedPoints polygon points, ordered
     */
    void sort_and_remove( Matrix< DDRMat > &aIntersectedPoints );

    /**
     * Computes the cross product of three points in 2D space, which is used to determine the orientation of the points.
     * @param[ in ] p1 First point
     * @param[ in ] p2 Second point
     * @param[ in ] p3 Third point
     * @return The cross product value, which indicates the orientation of the points.
     *         - If the value is positive, the points are oriented counter-clockwise.
     *         - If the value is negative, the points are oriented clockwise.
     *         - If the value is zero, the points are collinear.
     */
    real cross_tri( const Matrix< DDRMat > &p1, const Matrix< DDRMat > &p2, const Matrix< DDRMat > &p3 );

    /**
     * Computes the orientation of three points in 2D space.
     * @param[ in ] p First point
     * @param[ in ] q Second point
     * @param[ in ] r Third point
     * @return The orientation value, which indicates the relative position of the points.
     *         - If the value is positive, the points are oriented counter-clockwise.
     *         - If the value is negative, the points are oriented clockwise.
     *         - If the value is zero, the points are collinear.
     */
    real orientation( const Matrix< DDRMat > &p, const Matrix< DDRMat > &q, const Matrix< DDRMat > &r );

    /**
     * Checks if two line segments (edges) intersect properly.
     * @param[ in ] p1 First endpoint of the first line segment
     * @param[ in ] q1 Second endpoint of the first line segment
     * @param[ in ] p2 First endpoint of the second line segment
     * @param[ in ] q2 Second endpoint of the second line segment
     * @return true if the line segments intersect properly, false otherwise.
     */
    bool proper_edge_intersect( const Matrix< DDRMat > &p1, const Matrix< DDRMat > &q1, const Matrix< DDRMat > &p2, const Matrix< DDRMat > &q2 );

    /**
     * Checks if a point lies strictly inside a triangle in 2D space.
     * @param[ in ] aTriangle Coordinates of the triangle (3 vertices)
     * @param[ in ] aPoint Coordinates of the point to check
     * @return true if the point lies strictly inside the triangle, false otherwise.
     * @note The function uses the orientation of the triangle's edges and the point to determine if the point is inside the triangle.
     */
    bool strictly_contains( const Matrix< DDRMat > &aTriangle, const Matrix< DDRMat > &aPoint );

    /**
     * Checks if two triangles strictly overlap in 2D space.
     * @param[ in ] aFirstTriangle Coordinates of the first triangle (3 vertices)
     * @param[ in ] aSecondTriangle Coordinates of the second triangle (3 vertices)
     * @return true if the triangles strictly overlap (i.e., they have a non-empty intersection), false otherwise.
     * @note The function checks for edge-edge intersections and containment of one triangle's vertices within the other triangle.
     */
    bool triangles_strictly_overlap(
            const Matrix< DDRMat > &aFirstTriangle,
            const Matrix< DDRMat > &aSecondTriangle );
}    // namespace moris::mig

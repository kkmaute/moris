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


    // TODO @bc DOCUMENTATION AND VARIABLE RENAMING FOR ALL

    real cross_tri( const Matrix< DDRMat > &p1, const Matrix< DDRMat > &p2, const Matrix< DDRMat > &p3 );

    real orientation( const Matrix< DDRMat > &p, const Matrix< DDRMat > &q, const Matrix< DDRMat > &r );

    bool proper_edge_intersect( const Matrix< DDRMat > &p1, const Matrix< DDRMat > &q1, const Matrix< DDRMat > &p2, const Matrix< DDRMat > &q2 );

    bool strictly_contains( const Matrix< DDRMat > &aTriangle, const Matrix< DDRMat > &aPoint );

    bool triangles_strictly_overlap(
            const Matrix< DDRMat > &aFirstTriangle,
            const Matrix< DDRMat > &aSecondTriangle );
}    // namespace moris::mig

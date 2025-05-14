/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_MIG_Triangle_Intersect.cpp
 *
 */

#include "fn_MIG_Triangle_Intersect.hpp"

#include "fn_join_horiz.hpp"
#include "fn_inv.hpp"
#include "fn_sort.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"


namespace moris::mig
{
    void
    Intersect(
            Matrix< DDRMat > const & aFirstTRICoords,
            Matrix< DDRMat > const & aSecondTRICoords,
            Matrix< DDRMat >&        aIntersectedPoints )
    {
        edge_intersect( aFirstTRICoords, aSecondTRICoords, aIntersectedPoints );

        Matrix< DDRMat > P1;

        find_vertices_inside_triangle( aFirstTRICoords, aSecondTRICoords, P1 );

        aIntersectedPoints = join_horiz( aIntersectedPoints, P1 );

        find_vertices_inside_triangle( aSecondTRICoords, aFirstTRICoords, P1 );

        aIntersectedPoints = join_horiz( aIntersectedPoints, P1 );

        sort_and_remove( aIntersectedPoints );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool triangles_intersect(
            Matrix< DDRMat > const & aFirstTRICoords,
            Matrix< DDRMat > const & aSecondTRICoords )
    {
        Matrix< DDRMat > tIntersectionPoints;

        // check if there are any intersections along any edge of the triangles
        edge_intersect( aFirstTRICoords, aSecondTRICoords, tIntersectionPoints );

        if ( tIntersectionPoints.n_cols() > 0 )
        {
            return true;
        }

        // Check if the first triangle has any vertices inside the second triangle
        find_vertices_inside_triangle( aFirstTRICoords, aSecondTRICoords, tIntersectionPoints );
        if ( tIntersectionPoints.n_cols() > 0 )
        {
            return true;
        }

        // Check if the second triangle has any vertices inside the first triangle
        find_vertices_inside_triangle( aSecondTRICoords, aFirstTRICoords, tIntersectionPoints );
        if ( tIntersectionPoints.n_cols() > 0 )
        {
            return true;
        }
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    edge_intersect(
            Matrix< DDRMat > const & aFirstTRICoords,
            Matrix< DDRMat > const & aSecondTRICoords,
            Matrix< DDRMat >&        aIntersectedPoints )
    {
        uint tNumIntersections = 0;

        for ( uint i = 0; i < 3; i++ )
        {
            for ( uint j = 0; j < 3; j++ )
            {
                // form the matrix on RHS
                Matrix< DDRMat > b = aSecondTRICoords.get_column( j ) - aFirstTRICoords.get_column( i );

                // matrix on the LHS
                Matrix< DDRMat > A( 2, 2 );
                A.set_column( 0, aFirstTRICoords.get_column( ( i + 1 ) % 3 ) - aFirstTRICoords.get_column( i ) );
                A.set_column( 1, -aSecondTRICoords.get_column( ( j + 1 ) % 3 ) + aSecondTRICoords.get_column( j ) );

                // solve the system
                if ( std::abs( A( 0, 0 ) * A( 1, 1 ) - A( 0, 1 ) * A( 1, 0 ) ) > 0.0000001 )
                {
                    Matrix< DDRMat > R = inv( A ) * b;

                    moris::real eps = 0.01;

                    // Intersection Condition
                    if ( R( 0 ) >= 0 and ( R( 0 ) - 1 ) <= eps and R( 1 ) >= 0 and ( R( 1 ) - 1 ) <= eps )
                    {
                        // grow the matrix to insert the new intersection
                        aIntersectedPoints.resize( 2, tNumIntersections + 1 );

                        aIntersectedPoints.get_column( tNumIntersections ) = aFirstTRICoords.get_column( i ) + R( 0 ) * A.get_column( 0 );

                        // increase number of intersection by 1
                        tNumIntersections++;
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void find_vertices_inside_triangle(
            Matrix< DDRMat > const & aFirstTRICoords,
            Matrix< DDRMat > const & aSecondTRICoords,
            Matrix< DDRMat >&        aIntersectedPoints )
    {
        uint tNumIntersections = 0;

        // Interior points
        Matrix< moris::DDRMat > v0 = aSecondTRICoords.get_column( 1 ) - aSecondTRICoords.get_column( 0 );
        Matrix< moris::DDRMat > v1 = aSecondTRICoords.get_column( 2 ) - aSecondTRICoords.get_column( 0 );

        // Baricenteric Coordinates
        real d00 = dot( v0, v0 );
        real d01 = dot( v0, v1 );
        real d11 = dot( v1, v1 );

        real id = 1.0 / ( d00 * d11 - d01 * d01 );

        for ( uint i = 0; i < 3; i++ )
        {
            Matrix< moris::DDRMat > v2 = aFirstTRICoords.get_column( i ) - aSecondTRICoords.get_column( 0 );

            real d02 = dot( v0, v2 );
            real d12 = dot( v1, v2 );

            real u = ( d11 * d02 - d01 * d12 ) * id;
            real v = ( d00 * d12 - d01 * d02 ) * id;

            real eps = 0.001;

            if ( u >= -eps and v >= -eps and ( u + v <= 1 + eps ) )
            {
                aIntersectedPoints.resize( 2, tNumIntersections + 1 );

                aIntersectedPoints.get_column( tNumIntersections ) = aFirstTRICoords.get_column( i );

                tNumIntersections++;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    sort_and_remove( Matrix< DDRMat >& aIntersectedPoints )
    {
        real eps = 0.0001;

        uint tNumColumn = aIntersectedPoints.n_cols();

        if ( tNumColumn > 0 )
        {
            // C
            Matrix< moris::DDRMat > C( 2, 1 );
            C.get_row( 0 ) = sum( aIntersectedPoints.get_row( 0 ) ) / tNumColumn;
            C.get_row( 1 ) = sum( aIntersectedPoints.get_row( 1 ) ) / tNumColumn;

            Matrix< moris::DDRMat > ao( 1, tNumColumn );

            // order polygon corners counter
            for ( uint i = 0; i < tNumColumn; i++ )
            {
                Matrix< moris::DDRMat > d = aIntersectedPoints.get_column( i ) - C;

                ao( i ) = std::atan2( d( 1 ), d( 0 ) );
            }

            Matrix< moris::DDRMat > aoSorted;

            Vector< uint > IdMatrix( tNumColumn );

            moris::sort( ao, aoSorted, "ascend", 1 );

            // find the index matrix
            for ( uint i = 0; i < tNumColumn; i++ )
            {
                for ( uint j = 0; j < tNumColumn; j++ )
                {
                    if ( ao( j ) == aoSorted( i ) )
                    {
                        IdMatrix( i ) = j;
                        break;
                    }
                }
            }

            // sort the points based on IdMatrix
            Matrix< moris::DDRMat > tSortedIntersectedPoints( 2, tNumColumn );

            for ( uint i = 0; i < tNumColumn; i++ )
            {
                tSortedIntersectedPoints.get_column( i ) = aIntersectedPoints.get_column( IdMatrix( i ) );
            }

            // remove duplicates of the point
            uint i = 0;
            uint j = 1;

            while ( j < tNumColumn )
            {
                if ( norm( tSortedIntersectedPoints.get_column( i ) - tSortedIntersectedPoints.get_column( j ) ) > eps )
                {
                    i++;

                    tSortedIntersectedPoints.get_column( i ) = tSortedIntersectedPoints.get_column( j );

                    j++;
                }
                else
                {
                    j++;
                }
            }

            tSortedIntersectedPoints.resize( 2, i + 1 );

            aIntersectedPoints = tSortedIntersectedPoints;
        }
    }

    real cross_tri( const Matrix< DDRMat >& p1, const Matrix< DDRMat >& p2, const Matrix< DDRMat >& p3 )
    {
        real x1 = p2( 0 ) - p1( 0 );
        real y1 = p2( 1 ) - p1( 1 );
        real x2 = p3( 0 ) - p1( 0 );
        real y2 = p3( 1 ) - p1( 1 );
        return x1 * y2 - y1 * x2;
    }

    // Orientation of triplet (p, q, r)
    real orientation( const Matrix< DDRMat >& p, const Matrix< DDRMat >& q, const Matrix< DDRMat >& r )
    {
        return cross_tri( p, q, r );
    }

    // Check if two edges properly intersect (excluding touching at endpoints)
    bool proper_edge_intersect( const Matrix< DDRMat >& p1, const Matrix< DDRMat >& q1, const Matrix< DDRMat >& p2, const Matrix< DDRMat >& q2 )
    {
        real o1 = orientation( p1, q1, p2 );
        real o2 = orientation( p1, q1, q2 );
        real o3 = orientation( p2, q2, p1 );
        real o4 = orientation( p2, q2, q1 );

        return ( o1 * o2 < 0 ) && ( o3 * o4 < 0 );
    }

    // Check if point is strictly inside triangle (not on edge)
    bool strictly_contains( const Matrix< DDRMat >& aTriangle, const Matrix< DDRMat >& aPoint )
    {
        real tOrient1 = orientation( aTriangle.get_column( 0 ), aTriangle.get_column( 1 ), aPoint );
        real tOrient2 = orientation( aTriangle.get_column( 1 ), aTriangle.get_column( 2 ), aPoint );
        real tOrient3 = orientation( aTriangle.get_column( 2 ), aTriangle.get_column( 0 ), aPoint );
        return ( tOrient1 > 0 and tOrient2 > 0 and tOrient3 > 0 ) or ( tOrient1 < 0 and tOrient2 < 0 and tOrient3 < 0 );
    }

    // Check if triangles strictly overlap (excluding touching)
    bool triangles_strictly_overlap(
            const Matrix< DDRMat >& aFirstTriangle,
            const Matrix< DDRMat >& aSecondTriangle )
    {
        MORIS_ASSERT( aFirstTriangle.n_cols() == 3, "First triangle must have 3 vertices." );
        MORIS_ASSERT( aSecondTriangle.n_cols() == 3, "Second triangle must have 3 vertices." );

        // Edge-edge intersection
        for ( uint iFirstTriPoint = 0; iFirstTriPoint < 3; ++iFirstTriPoint )
        {
            const Matrix< DDRMat >& p1 = aFirstTriangle.get_column( iFirstTriPoint );
            const Matrix< DDRMat >& q1 = aFirstTriangle.get_column( ( iFirstTriPoint + 1 ) % 3 );
            for ( int iSecondTriPoint = 0; iSecondTriPoint < 3; ++iSecondTriPoint )
            {
                const Matrix< DDRMat >& p2 = aSecondTriangle.get_column( iSecondTriPoint );
                const Matrix< DDRMat >& q2 = aSecondTriangle.get_column( ( iSecondTriPoint + 1 ) % 3 );
                if ( proper_edge_intersect( p1, q1, p2, q2 ) )
                {
                    return true;
                }
            }
        }

        // Containment
        for ( uint iPoint = 0; iPoint < 3; ++iPoint )
        {
            const Matrix< DDRMat >& tPoint = aFirstTriangle.get_column( iPoint );
            if ( strictly_contains( aSecondTriangle, tPoint ) )
            {
                return true;
            }
        }
        for ( uint iPoint = 0; iPoint < 3; ++iPoint )
        {
            const Matrix< DDRMat >& tPoint = aSecondTriangle.get_column( iPoint );
            if ( strictly_contains( aFirstTriangle, tPoint ) )
            {
                return true;
            }
        }

        return false;
    }
}    // namespace moris::mig
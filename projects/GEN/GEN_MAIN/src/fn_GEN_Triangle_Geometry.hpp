/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_Triangle_Geometry.hpp
 *
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMENG_FN_GEN_TRIANGLE_GEOMETRY_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMENG_FN_GEN_TRIANGLE_GEOMETRY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_comp_abs.hpp"
#include "fn_tet_volume.hpp"
#include "fn_tri_area.hpp"

namespace moris
{
    namespace ge
    {
        /*
         * Compute the surface normal of a triangle (Note this hasn't been fully optimized)
         * @param[in] aTriangleNodes - Nodes of the triangle
         */
        inline void
        compute_tri_surface_normal(
                moris::Matrix< moris::IndexMat > const & aTriangleNodes,
                moris::Matrix< moris::DDRMat > const &   aNodeCoordinates,
                moris::Matrix< moris::DDRMat >&          aSurfaceNormal,
                bool                                     aPosDirection = false )
        {
            moris::Matrix< moris::DDRMat > tVec1( 1, 3, 0 );
            moris::Matrix< moris::DDRMat > tVec2( 1, 3, 0 );

            moris::moris_index tN1 = aTriangleNodes( 0, 0 );
            moris::moris_index tN2 = aTriangleNodes( 0, 1 );
            moris::moris_index tN3 = aTriangleNodes( 0, 2 );

            tVec1( 0, 0 ) = aNodeCoordinates( tN3, 0 ) - aNodeCoordinates( tN1, 0 );
            tVec1( 0, 1 ) = aNodeCoordinates( tN3, 1 ) - aNodeCoordinates( tN1, 1 );
            tVec1( 0, 2 ) = aNodeCoordinates( tN3, 2 ) - aNodeCoordinates( tN1, 2 );

            tVec2( 0, 0 ) = aNodeCoordinates( tN2, 0 ) - aNodeCoordinates( tN1, 0 );
            tVec2( 0, 1 ) = aNodeCoordinates( tN2, 1 ) - aNodeCoordinates( tN1, 1 );
            tVec2( 0, 2 ) = aNodeCoordinates( tN2, 2 ) - aNodeCoordinates( tN1, 2 );

            // Compute cross products
            aSurfaceNormal( 0, 0 ) = tVec1( 0, 1 ) * tVec2( 0, 2 ) - tVec1( 0, 2 ) * tVec2( 0, 1 );
            aSurfaceNormal( 1, 0 ) = tVec1( 0, 2 ) * tVec2( 0, 0 ) - tVec1( 0, 0 ) * tVec2( 0, 2 );
            aSurfaceNormal( 2, 0 ) = tVec1( 0, 0 ) * tVec2( 0, 1 ) - tVec1( 0, 1 ) * tVec2( 0, 0 );

            // Normalized
            moris::real tLenSquared3 = std::pow( aSurfaceNormal( 0, 0 ), 2 ) + std::pow( aSurfaceNormal( 1, 0 ), 2 ) + std::pow( aSurfaceNormal( 2, 0 ), 2 );
            moris::real tLen3        = std::pow( tLenSquared3, 0.5 );

            MORIS_ASSERT( tLen3 > MORIS_REAL_MIN, "Dividing by near zero value" );

            moris::real tLenInv = 1.0 / tLen3;

            if ( !aPosDirection )
            {
                aSurfaceNormal = tLenInv * aSurfaceNormal;
            }
            else
            {
                aSurfaceNormal = tLenInv * comp_abs( aSurfaceNormal );
            }
        }

        //-----------------------------------------------------------------------

        inline moris::real
        compute_volume_for_multiple_tets(
                moris::Matrix< moris::DDRMat > const &   aAllNodeCoords,
                moris::Matrix< moris::IndexMat > const & aElementToNodeConnectivity )
        {
            size_t tNumTets = aElementToNodeConnectivity.n_rows();

            moris::Matrix< moris::DDRMat > tCoords( aElementToNodeConnectivity.n_cols(), 3 );
            moris::Matrix< moris::DDRMat > tCoordRow( 1, 3 );

            moris::real tTotalVolume = 0;
            moris::real tSingleVolume;
            size_t      k = 0;

            for ( size_t i = 0; i < tNumTets; i++ )
            {
                for ( size_t j = 0; j < aElementToNodeConnectivity.n_cols(); j++ )
                {
                    tCoordRow = aAllNodeCoords.get_row( aElementToNodeConnectivity( i, j ) );
                    tCoords.set_row( k, tCoordRow );
                    k++;
                }
                k             = 0;
                tSingleVolume = xtk::vol_tetrahedron( tCoords );
                //        std::cout<<"    i = "<<i<<" tSingleVolume = "<<tSingleVolume<<std::endl;
                tTotalVolume += tSingleVolume;
            }

            return tTotalVolume;
        }

        //-----------------------------------------------------------------------

        template< typename Real_Matrix >
        inline moris::Matrix< Real_Matrix >
        compute_val_at_tet_centroid( moris::Matrix< Real_Matrix > const & aNodeVals )
        {
            Real_Matrix tVal = 0.25 * ( aNodeVals );

            moris::Matrix< Real_Matrix > tValMat( tVal );

            return tValMat;
        }
    }    // namespace ge

    //-----------------------------------------------------------------------

    inline moris::real
    compute_area_for_multiple_triangles(
            moris::Matrix< moris::DDRMat > const &   aAllNodeCoords,
            moris::Matrix< moris::IndexMat > const & aElementToNodeConnectivity )
    {
        size_t tNumTriangles     = aElementToNodeConnectivity.n_rows();
        size_t tSpacialDimension = aElementToNodeConnectivity.n_cols();

        moris::Matrix< moris::DDRMat > tCoords( aElementToNodeConnectivity.n_cols(), 2 );
        moris::Matrix< moris::DDRMat > tCoordRow( 1, 2 );

        moris::real tTotalArea = 0;
        moris::real tSingleArea;
        size_t      k = 0;

        for ( size_t i = 0; i < tNumTriangles; i++ )
        {
            for ( size_t j = 0; j < tSpacialDimension; j++ )
            {
                tCoordRow = aAllNodeCoords.get_row( aElementToNodeConnectivity( i, j ) );
                tCoords.set_row( k, tCoordRow );
                k++;
            }
            k           = 0;
            tSingleArea = xtk::area_tri_2D( tCoords );
            tTotalArea += tSingleArea;
        }

        return tTotalArea;
    }
}    // namespace moris

#endif /* PROJECTS_GEN_SRC_NEW_GEOMENG_FN_GEN_TRIANGLE_GEOMETRY_HPP_ */

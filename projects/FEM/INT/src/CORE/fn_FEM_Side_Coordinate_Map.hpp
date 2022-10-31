/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Side_Coordinate_Map.hpp
 *
 */

#ifndef SRC_FEM_FN_FEM_SIDE_COORDINATE_MAP_HPP_
#define SRC_FEM_FN_FEM_SIDE_COORDINATE_MAP_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    //------------------------------------------------------------------------------
    namespace fem
    {

        //------------------------------------------------------------------------------
        // function that maps local coordinates on sides of hex (side is quad) and
        // tets (side is tri) from master to slave side
        //
        // function overwrites only spatial coordinates assuming that spatial coordinates
        // are first entries of vector

        inline void
        side_coordinate_map(
                const mtk::Geometry_Type       aGeometryType,
                const moris_index              aSlaveNode,
                const moris::Matrix< DDRMat >& aMasterCoordinates,
                moris::Matrix< DDRMat >&       aSlaveCoordinates )
        {
            // switch on geometry type
            switch ( aGeometryType )
            {
                case mtk::Geometry_Type::QUAD:
                {
                    MORIS_ASSERT( aMasterCoordinates.n_rows() > 1,
                            "side_coordinate_map - master coordinate vector has insufficient number of rows." );

                    moris::Matrix< DDRMat > tRotationMatrix;

                    switch ( aSlaveNode )
                    {
                        case 0:
                            tRotationMatrix = { { 0.0, 1.0 }, { 1.0, 0.0 } };
                            break;
                        case 1:
                            tRotationMatrix = { { -1.0, 0.0 }, { 0.0, 1.0 } };
                            break;
                        case 2:
                            tRotationMatrix = { { 0.0, -1.0 }, { -1.0, 0.0 } };
                            break;
                        case 3:
                            tRotationMatrix = { { 1.0, 0.0 }, { 0.0, -1.0 } };
                            break;
                        default:
                            MORIS_ERROR( false, " side_coordinate_map - unknown slave node " );
                    }

                    // write mapped coordinates onto first two rows
                    aSlaveCoordinates( { 0, 1 } ) = tRotationMatrix * aMasterCoordinates( { 0, 1 } );

                    break;
                }

                case mtk::Geometry_Type::TRI:
                {
                    MORIS_ASSERT( aMasterCoordinates.n_rows() > 1,
                            "side_coordinate_map - master coordinate vector has insufficient number of rows." );

                    moris::Matrix< DDRMat > tRotationMatrix;
                    moris::Matrix< DDRMat > tOffsetVecor;

                    switch ( aSlaveNode )
                    {
                        case 0:
                            tRotationMatrix = { { 1.0, 0.0 }, { -1.0, -1.0 } };
                            tOffsetVecor    = { { 0.0 }, { 1.0 } };
                            break;
                        case 1:
                            tRotationMatrix = { { 0.0, 1.0 }, { 1.0, 0.0 } };
                            tOffsetVecor    = { { 0.0 }, { 0.0 } };
                            break;
                        case 2:
                            tRotationMatrix = { { -1.0, -1.0 }, { 0.0, 1.0 } };
                            tOffsetVecor    = { { 1.0 }, { 0.0 } };
                            break;
                        default:
                            MORIS_ASSERT( false, " side_coordinate_map - unknown slave node " );
                    }

                    // write mapped coordinates onto first two rows
                    aSlaveCoordinates( { 0, 1 } ) = tOffsetVecor + tRotationMatrix * aMasterCoordinates( { 0, 1 } );

                    break;
                }

                case mtk::Geometry_Type::LINE:
                {
                    // write mapped coordinate onto first row
                    aSlaveCoordinates( 0 ) = -aMasterCoordinates( 0 );
                    ;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, " side_coordinate_map - unknown geometry type " );
                }
            }
        }

        //------------------------------------------------------------------------------

        inline moris::Matrix< DDRMat >
        side_coordinate_map(
                const mtk::Geometry_Type       aGeometryType,
                const moris_index              aSlaveNode,
                const moris::Matrix< DDRMat >& aMasterCoordinates )
        {
            // check for proper dimensions of master coordinate vector
            MORIS_ASSERT( aMasterCoordinates.n_cols() == 1,
                    "side_coordinate_map - master coordinate vector is not column vector." );

            // initialize slave coordinate vector with master coordinate vector,
            moris::Matrix< DDRMat > tSlaveCoordinates = aMasterCoordinates;

            // evaluate map
            side_coordinate_map( aGeometryType, aSlaveNode, aMasterCoordinates, tSlaveCoordinates );

            // return slave coordinate vector
            return tSlaveCoordinates;
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_FN_FEM_SIDE_COORDINATE_MAP_HPP_ */

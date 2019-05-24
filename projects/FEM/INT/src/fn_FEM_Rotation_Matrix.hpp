/*
 * cl_FEM_fn_fem_rotation_matrix.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_FN_FEM_ROTATION_MATRIX_HPP_
#define SRC_FEM_FN_FEM_ROTATION_MATRIX_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
//------------------------------------------------------------------------------
    namespace fem
    {
//------------------------------------------------------------------------------

//        Matrix< DDRMat> rotation_matrix( mtk::Geometry_Type aGeometryType,
//                                         moris_index        aMasterSideOrdinal,
//                                         moris_index        aSlaveSideOrdinal,
//                                         moris_index        aSlaveNode )
        Matrix< DDRMat> rotation_matrix( mtk::Geometry_Type aGeometryType,
                                         moris_index        aSlaveNode )
        {
            // init the rotation matrix
            Matrix < DDRMat > tR;

            // switch on geometry type
            switch( aGeometryType )
            {
                case( mtk::Geometry_Type::HEX ):
                {
                    switch( aSlaveNode )
                    {
                        case( 0 ):
                            tR = {{ 0.0, 1.0 },{ 1.0, 0.0 }};
                            break;
                        case( 1 ):
                            tR = {{ -1.0, 0.0 },{ 0.0, 1.0 }};
                            break;
                        case( 2 ):
                            tR = {{ 0.0, -1.0 },{ -1.0, 0.0 }};
                            break;
                        case( 3 ):
                            tR = {{ 1.0, 0.0 },{ 0.0, -1.0 }};
                            break;
                        default:
                            MORIS_ERROR( false, " rotation_matrix - unknown slave node ");
                            break;
                    }
                    break;
                }
                case( mtk::Geometry_Type::TET ):
                {
                    switch( aSlaveNode )
                    {
                        case( 0 ):
                            tR = {{ 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0 }};
                            break;
                        case( 1 ):
                            tR = {{ 0.0, 1.0, 0.0 }, {1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }};
                            break;
                        case( 2 ):
                            tR = {{ 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }};
                            break;
                        default:
                            MORIS_ASSERT( false, " rotation_matrix - unknown slave node " );
                            break;
                    }
                    break;
                }
                case( mtk::Geometry_Type::QUAD ):
                case( mtk::Geometry_Type::TRI ):
                {
                    tR = {{ -1.0 }};
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, " rotation_matrix - unknown geometry type ");
                    break;
                }
            }
            return tR;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_FN_FEM_ROTATION_MATRIX_HPP_ */

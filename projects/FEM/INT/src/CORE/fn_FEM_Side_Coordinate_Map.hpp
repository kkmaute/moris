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

        void rotation_matrix( 
                mtk::Geometry_Type        aGeometryType,
                moris_index               aSlaveNode,
                moris::Matrix< DDRMat > & aRotationMatrix )
        {
            // switch on geometry type
            switch( aGeometryType )
            {
                case mtk::Geometry_Type::QUAD :
                {
                    switch( aSlaveNode )
                    {
                        case 0 :
                            aRotationMatrix = {{ 0.0, 1.0 },{ 1.0, 0.0 }};
                            break;
                        case 1 :
                            aRotationMatrix = {{ -1.0, 0.0 },{ 0.0, 1.0 }};
                            break;
                        case 2 :
                            aRotationMatrix = {{ 0.0, -1.0 },{ -1.0, 0.0 }};
                            break;
                        case 3 :
                            aRotationMatrix = {{ 1.0, 0.0 },{ 0.0, -1.0 }};
                            break;
                        default:
                            MORIS_ERROR( false, " rotation_matrix - unknown slave node ");
                    }
                    break;
                }

                case mtk::Geometry_Type::TRI :
                {
                    switch( aSlaveNode )
                    {
                        case 0 :
                            aRotationMatrix = {{ 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0 }};
                            break;
                        case 1 :
                            aRotationMatrix = {{ 0.0, 1.0, 0.0 }, {1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }};
                            break;
                        case 2 :
                            aRotationMatrix = {{ 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }};
                            break;
                        default:
                            MORIS_ASSERT( false, " rotation_matrix - unknown slave node " );
                    }
                    break;
                }

                case mtk::Geometry_Type::LINE :
                {
                    aRotationMatrix = {{ -1.0 }};
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, " rotation_matrix - unknown geometry type ");
                }
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_FN_FEM_ROTATION_MATRIX_HPP_ */

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_IWG_Compressible_NS_Tools.cpp
 *
 */

#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Matrix< DDRMat >
        unfold_flat_tensor( const Matrix< DDRMat >& aFlattenedTensor )
        {
            // get length of the flattened tensor
            uint tLength = aFlattenedTensor.length();

            // clang-format off
            switch ( tLength )
            {
                // 2D: convert 3x1 flattened tensor back to 2x2 matrix
                case 3 :
                {
                    return {
                            { aFlattenedTensor( 0 ), aFlattenedTensor( 2 ) },
                            { aFlattenedTensor( 2 ), aFlattenedTensor( 1 ) } };
                    break;
                }

                // 3D: convert 6x1 flattened tensor back to 3x3 matrix
                case 6 :
                {
                    return {
                            { aFlattenedTensor( 0 ), aFlattenedTensor( 5 ), aFlattenedTensor( 4 ) },
                            { aFlattenedTensor( 5 ), aFlattenedTensor( 1 ), aFlattenedTensor( 3 ) },
                            { aFlattenedTensor( 4 ), aFlattenedTensor( 3 ), aFlattenedTensor( 2 ) } };
                    break;
                }

                default:
                {
                    MORIS_ERROR( false,
                            "fn_FEM_IWG_Compressible_NS::unfold_flat_tensor - aFlattenedTensor must be a vector of length 3 (2D) or 6 (3D)." );
                    return { { 0.0 } };
                    break;
                }
            }
            // clang-format on
        }

        //------------------------------------------------------------------------------

        uint
        convert_index_pair_to_flat(
                const uint aI,
                const uint aJ,
                const uint aNumSpaceDims )
        {
            // check that indices are not out of bounds
            MORIS_ASSERT( ( aI < aNumSpaceDims ) and ( aJ < aNumSpaceDims ),
                    "fn_FEM_IWG_Compressible_NS::convert_index_pair_to_flat() - indices out of bounds." );

            // initialize return value
            uint tFlatIndex = 0;

            // clang-format off
            // two cases: 2D and 3D
            switch ( aNumSpaceDims )
            {

                // 2D
                case 2 :
                {
                    if ( aI + aJ == 0 )
                    {
                        tFlatIndex = 0;
                    }
                    else if ( aI + aJ == 2 )
                    {
                        tFlatIndex = 1;
                    }
                    else
                    {
                        tFlatIndex = 2;
                    }
                    break;
                }

                // 3D
                case 3 :
                {
                    if ( aI + aJ == 0 )
                    {
                        tFlatIndex = 0;
                    }
                    else if ( aI + aJ == 2 )
                    {
                        tFlatIndex = 1;
                    }
                    else if ( aI + aJ == 5 )
                    {
                        tFlatIndex = 3;
                    }
                    else if ( aI + aJ == 3 )
                    {
                        tFlatIndex = 5;
                    }
                    else if ( aI == 2 )
                    {
                        tFlatIndex = 2;
                    }
                    else
                    {
                        tFlatIndex = 4;
                    }
                    break;
                }

                default :
                {
                    MORIS_ERROR( false,
                            "fn_FEM_IWG_Compressible_NS::convert_index_pair_to_flat() - "
                            "number of spatial dimensions not supported (only 2D or 3D supported)." );
                    break;
                }

            } // end switch statement
            // clang-format on

            // return index
            return tFlatIndex;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * fn_XTK_match_normal_to_side_ordinal.hpp  
 * 
 */
#ifndef SRC_fn_XTK_match_normal_to_side_ordinal
#define SRC_fn_XTK_match_normal_to_side_ordinal

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "typedefs.hpp"

using namespace moris;

namespace xtk
{

    //-------------------------------------------------------------------------------------

    /**
     * @brief checks if a normal corresponds to one of the global axis directions and returns the corresponding ordinal for a QUAD/HEX
     * NOTE: instead of failing, MORIS_INDEX_MAX is returned if the normal does not match any coordinate direction.
     * 
     * @param aNormal normal for which the ordinal is to be found
     * @return moris_index ordinal for a QUAD/HEX the coordinate direction corresponds to
     */
    moris_index 
    match_normal_to_side_ordinal( Matrix< DDRMat > const & aNormal )
    {
        // define a tolerance 
        real tEps = 1.0e-6;

        // normalize the normal
        Matrix< DDRMat > tNormal = aNormal / norm( aNormal );

        // convert normal to row vector if it isn't already
        if ( tNormal.n_cols() == 1 )
        {
            tNormal = trans( tNormal );
        }
        
        // case: 2D
        if ( aNormal.numel() == 2 )
        {
            const Cell< Matrix< DDRMat > > tOrdinalNormals = { 
                    { {  0.0, -1.0 } },   // 0 
                    { {  1.0,  0.0 } },   // 1
                    { {  0.0,  1.0 } },   // 2
                    { { -1.0,  0.0 } } }; // 3

            // check if any of the ordinals match
            for( uint iOrdinal = 0; iOrdinal < 4; iOrdinal++ )
            {
                real tNorm = norm( tNormal - tOrdinalNormals( iOrdinal ) );
                if ( tNorm < tEps )
                {
                    return iOrdinal;
                }
            }

            // if non of the ordinals match, return a false
            return MORIS_INDEX_MAX;
        }

        // case: 3D
        else if ( aNormal.numel() == 3 ) 
        {
            const Cell< Matrix< DDRMat > > tOrdinalNormals = { 
                    { {  0.0, -1.0,  0.0 } },   // 0 
                    { {  1.0,  0.0,  0.0 } },   // 1
                    { {  0.0,  1.0,  0.0 } },   // 2
                    { { -1.0,  0.0,  0.0 } },   // 3
                    { {  0.0,  0.0, -1.0 } },   // 4
                    { {  0.0,  0.0,  1.0 } } }; // 5

            // check if any of the ordinals match
            for( uint iOrdinal = 0; iOrdinal < 6; iOrdinal++ )
            {
                real tNorm = norm( tNormal - tOrdinalNormals( iOrdinal ) );
                if ( tNorm < tEps )
                {
                    return iOrdinal;
                }
            }

            // if non of the ordinals match, return a false
            return MORIS_INDEX_MAX;

        }

        // dimension unknown
        else
        {
            MORIS_ERROR( false, "xtk::match_normal_to_side_ordinal() - Normal provided is neither 2D nor 3D." );
            return MORIS_INDEX_MAX;
        }
    }

    //-------------------------------------------------------------------------------------

}    // namespace xtk


#endif /* fn_XTK_match_normal_to_side_ordinal.hpp */

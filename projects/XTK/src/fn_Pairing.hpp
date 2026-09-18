/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Pairing.hpp
 *
 */

#ifndef SRC_TOOLS_FN_PAIRING_HPP_
#define SRC_TOOLS_FN_PAIRING_HPP_

#include "cl_Matrix.hpp"

namespace moris::xtk
{
    template< typename Integer >
    inline static Integer
    cantor_pairing(
            Integer const aPt1,
            Integer const aPt2 )
    {
        enum PairingType : uint
        {
            CANTOR,
            SZUDZIK
        };

        const PairingType tPairing = PairingType::CANTOR;

        Integer tID;

        switch ( tPairing )
        {
            case PairingType::CANTOR:
            {
                // Cantor pairing
                tID = ( ( aPt1 * aPt1 + 3 * aPt1 + 2 * aPt1 * aPt2 + aPt2 + aPt2 * aPt2 ) / 2 );
                break;
            }
            case PairingType::SZUDZIK:
            {
                // Szudzik Pairing
                if ( aPt1 > aPt2 )
                {
                    tID = (Integer)aPt1 * aPt1 + aPt1 + aPt2;
                }
                else
                {
                    tID = (Integer)aPt2 * aPt2 + aPt1;
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "cantor_pairing: Unknown pairing type" );
                break;
            }
        }

        return tID;
    }
}    // namespace moris::xtk

#endif /* SRC_TOOLS_FN_PAIRING_HPP_ */

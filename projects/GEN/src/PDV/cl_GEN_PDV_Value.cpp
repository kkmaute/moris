/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV_Value.cpp
 *
 */

#include "cl_GEN_PDV_Value.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    PDV_Value::PDV_Value( real aValue )
            : mValue( aValue )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool PDV_Value::is_active()
    {
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    real PDV_Value::get_value( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates )
    {
        return mValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > PDV_Value::get_sensitivities( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates )
    {
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > PDV_Value::get_determining_adv_ids( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates )
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

}

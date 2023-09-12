/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Pdv_Value.cpp
 *
 */

#include "cl_GEN_Pdv_Value.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Value::Pdv_Value(real aValue)
        : mValue(aValue)
        {
            mIsActive = false;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Pdv_Value::get_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return mValue;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Pdv_Value::get_sensitivities(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return {{}};
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Pdv_Value::get_determining_adv_ids(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return {{}};
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


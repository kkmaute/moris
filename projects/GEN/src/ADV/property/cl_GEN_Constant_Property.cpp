/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Constant_Property.cpp
 *
 */

#include "cl_GEN_Constant_Property.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    real Constant_Property::get_field_value(const Matrix<DDRMat>& aCoordinates)
    {
        return *mVariables(0);
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>& Constant_Property::get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates)
    {
        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Constant_Property::get_dfield_dcoordinates(
            const Matrix<DDRMat>& aCoordinates,
            Matrix<DDRMat>&       aSensitivities)
    {
        MORIS_ERROR(false, "get_dfield_dcoordinates not implemented for constant property.");
    }

    //--------------------------------------------------------------------------------------------------------------
}


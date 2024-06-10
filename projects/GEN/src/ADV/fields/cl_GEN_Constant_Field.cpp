/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Constant_Field.cpp
 *
 */

#include "cl_GEN_Constant_Field.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Constant_Field::Constant_Field(
            const ADV&  aConstant,
            std::string aName )
            : Field_Analytic< 0 >( { aConstant }, std::move( aName ) )
    {
        mSensitivities = { { 1.0 } };
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Constant_Field::get_field_value( const Matrix< DDRMat >& aCoordinates )
    {
        return mADVHandler.get_variable( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Constant_Field::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
    {
        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Constant_Field::get_dfield_dcoordinates(
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        MORIS_ERROR( false, "get_dfield_dcoordinates not implemented for constant field." );
    }

    //--------------------------------------------------------------------------------------------------------------
}

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Analytic.cpp
 *
 */

#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field_Analytic::Field_Analytic()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Field_Analytic::~Field_Analytic()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Field_Analytic::get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            return this->get_field_value( aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Analytic::get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            return this->get_dfield_dadvs( aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Field_Analytic::get_dfield_dcoordinates(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities )
        {
            this->get_dfield_dcoordinates( aCoordinates, aSensitivities );
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

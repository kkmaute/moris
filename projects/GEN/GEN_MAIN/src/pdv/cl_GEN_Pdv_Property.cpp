/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Pdv_Property.cpp
 *
 */

#include "cl_GEN_Pdv_Property.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------

        Pdv_Property::Pdv_Property( std::shared_ptr< Property > aPropertyPointer )
                : mProperty( aPropertyPointer )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Property::~Pdv_Property()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Pdv_Property::get_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            return mProperty->get_field_value( aNodeIndex, aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Pdv_Property::get_sensitivities(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            return mProperty->get_dfield_dadvs( aNodeIndex, aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDSMat >
        Pdv_Property::get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            return mProperty->get_determining_adv_ids( aNodeIndex, aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

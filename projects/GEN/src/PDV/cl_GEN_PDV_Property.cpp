/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV_Property.cpp
 *
 */

#include "cl_GEN_PDV_Property.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------

        PDV_Property::PDV_Property( std::shared_ptr< Property > aPropertyPointer )
                : mProperty( aPropertyPointer )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        PDV_Property::~PDV_Property()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        PDV_Property::get_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            return mProperty->get_field_value( aNodeIndex, aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        PDV_Property::get_sensitivities(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            return mProperty->get_dfield_dadvs( aNodeIndex, aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDSMat >
        PDV_Property::get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            return mProperty->get_determining_adv_ids( aNodeIndex, aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

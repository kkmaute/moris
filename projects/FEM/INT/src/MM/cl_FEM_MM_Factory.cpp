/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_MM_Factory.cpp
 *
 */

#include "assert.hpp"
//FEM/INT/src
#include "cl_FEM_MM_Factory.hpp"
#include "cl_FEM_MM_Perfect_Gas.hpp"
//#include "cl_FEM_MM_Van_Der_Waals_Fluid.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------
    std::shared_ptr< Material_Model > MM_Factory::create_MM( fem::Material_Type aMaterialType )
    {

        switch ( aMaterialType )
        {
            case Material_Type::PERFECT_GAS:
                return std::make_shared< MM_Perfect_Gas >();

                // case Material_Type::VAN_DER_WAALS_FLUID :
                //     return std::make_shared< MM_Van_Der_Waals_Fluid >();

                default:
                    MORIS_ERROR( false, " MM_Factory::create_MM - No material type specified. " );
                    return nullptr;
            }
        }
        //------------------------------------------------------------------------------
}    // namespace moris::fem

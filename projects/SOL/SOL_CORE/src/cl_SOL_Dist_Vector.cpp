/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Dist_Vector.cpp
 *
 */

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
    namespace sol
    {
        //--------------------------------------------------------------------------------------------------------------

        Dist_Vector::Dist_Vector(
                bool aManageMap )
                : mManageMap( aManageMap )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Dist_Vector::~Dist_Vector()
        {
        }

        //--------------------------------------------------------------------------------------------------------------
    }    // namespace sol
}    // namespace moris

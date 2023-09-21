/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_User_Defined_Property.cpp
 *
 */

#include "cl_GEN_User_Defined_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Property::User_Defined_Property(
                Matrix<DDRMat>            aConstants,
                Field_Function            aFieldFunction,
                Property_Parameters aParameters)
                : User_Defined_Field(aConstants, aFieldFunction, aParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


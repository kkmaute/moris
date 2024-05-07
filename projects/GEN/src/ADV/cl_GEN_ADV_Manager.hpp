/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_ADV_Manager.hpp
 *
 */

#pragma once

#include "cl_Design_Variable.hpp"
#include "cl_GEN_ADV.hpp"
#include "cl_Matrix.hpp"

namespace moris::gen
{
    class ADV_Manager
    {
      // TODO manage data privately
      public:
        Vector< real > mADVs;
        Vector< real > mLowerBounds;
        Vector< real > mUpperBounds;

        /**
         * Adds a design variable to the ADV manager, if it is active.
         *
         * @param aDesignVariable Design variable
         */
        ADV create_adv( const Design_Variable& aDesignVariable );
    };
}

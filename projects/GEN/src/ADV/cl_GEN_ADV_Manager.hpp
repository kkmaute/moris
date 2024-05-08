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
      private:
        Vector< char > mParameterIDs;
        bool mParameterIDsFinalized = false;

      // TODO manage data privately
      public:
        Vector< real > mADVs;
        Vector< real > mLowerBounds;
        Vector< real > mUpperBounds;

        /**
         * Registers the given unique parameter IDs in the ADV manager's vector of IDs.
         *
         * @param aParameterIDs
         */
        void register_parameter_ids( const Vector< char >& aParameterIDs );

        /**
         * Finalizes the parameter IDs such that ADVs and bounds can be reserved properly.
         */
        void finalize_parameter_ids();

        /**
         * Adds a design variable to the ADV manager, if it is active.
         *
         * @param aDesignVariable Design variable
         */
        ADV create_adv( const Design_Variable& aDesignVariable );
    };
}

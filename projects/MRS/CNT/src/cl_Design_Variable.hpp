/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_Design_Variable.hpp
*
*/

#include "moris_typedefs.hpp"

#pragma once

namespace moris
{
    /**
     * Design variable class, for ADV parameters
     */
    class Design_Variable
    {
      private:
        real mValue;
        real mLowerBound;
        real mUpperBound;
        bool mIsConstant;

      public:
        /**
         * Constructor for a constant variable.
         *
         * @param aConstantValue
         */
        Design_Variable( real aConstantValue );

        /**
         * Constructor for a varying design variable.
         *
         * @param aLowerBound Lower bound
         * @param aInitialValue Initial value
         * @param aUpperBound Upper value
         */
        Design_Variable( real aLowerBound, real aInitialValue, real aUpperBound );

        /**
         * Gets if this design variable is constant or not.
         *
         * @return Is constant
         */
        bool is_constant();

        /**
         * Gets the value of this design variable.
         *
         * @return Current value
         */
        real get_value();

        /**
         * Gets the lower bound of this design variable.
         *
         * @return Lower bound
         */
        real get_lower_bound();

        /**
         * Gets the upper bound of this design variable.
         *
         * @return Upper bound
         */
        real get_upper_bound();
    };
}

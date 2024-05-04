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
         * @param aConstantValue Constant value for this variable
         */
        Design_Variable( real aConstantValue = 0.0 ); // NOLINT implicit conversion is fine here

        /**
         * Constructor for a varying design variable.
         *
         * @param aLowerBound Variable lower bound
         * @param aInitialValue Initial value for this variable
         * @param aUpperBound Variable upper bound
         */
        Design_Variable( real aLowerBound, real aInitialValue, real aUpperBound );

        /**
         * Gets if this design variable is constant or not.
         *
         * @return Is constant
         */
        bool is_constant() const;

        /**
         * Gets the value of this design variable.
         *
         * @return Current value
         */
        real get_value() const;

        /**
         * Gets the lower bound of this design variable.
         *
         * @return Lower bound
         */
        real get_lower_bound() const;

        /**
         * Gets the upper bound of this design variable.
         *
         * @return Upper bound
         */
        real get_upper_bound() const;
    };

    /**
     * Less than operator for design variables, comparing both the value and bounds.
     *
     * @param aLeft Left design variable argument
     * @param aRight Right design variable argument
     * @return If left is less than right
     */
    bool operator <( const Design_Variable& aLeft, const Design_Variable& aRight );

    /**
     * Greater than operator for design variables, comparing both the value and bounds.
     *
     * @param aLeft Left design variable argument
     * @param aRight Right design variable argument
     * @return If left is greater than right
     */
    bool operator >( const Design_Variable& aLeft, const Design_Variable& aRight );
}

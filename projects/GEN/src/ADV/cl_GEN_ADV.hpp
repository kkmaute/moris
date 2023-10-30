/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_GEN_ADV.hpp
*
*/

#pragma once

#include "typedefs.hpp"
#include "cl_Matrix.hpp"

namespace moris::sol
{
    class Dist_Vector;
}

namespace moris::ge
{
    class ADV
    {
      private:
        real* mValue;
        sint mID;

      public:
        /**
         * Constructor given an ADV vector and an ID
         *
         * @param aADVs ADV vector
         * @param aID ID in the ADV vector
         */
        ADV( Matrix< DDRMat >& aADVs, sint aADVId );

        /**
         * Constructor given an ADV vector and an ID
         *
         * @param aADVs ADV vector
         * @param aID ID in the ADV vector
         */
        ADV( sol::Dist_Vector* aADVs, sint aADVId );

        /**
         * Constructor for a constant ADV, with a given value
         *
         * @param aValue
         */
        explicit ADV( real aValue );

        /**
         * Copy constructor such that constant ADVs create a new value instead of just copying the pointer.
         *
         * @param aCopyADV ADV to copy
         */
        ADV( const ADV& aCopyADV );

        /**
         * Destructor, ensures that a constant variable is cleaned up
         */
        ~ADV();

        /**
         * Gets the value of this ADV
         *
         * @return
         */
        real get_value();

        /**
         * Replaces the value of this ADV, if it is a constant value. Does not replace a value in the full ADV vector.
         *
         * @param aNewValue New constant value
         */
        void replace_constant( real aNewValue );
    };
}
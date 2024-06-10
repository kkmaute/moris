/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV_Value.hpp
 *
 */

#pragma once

#include "cl_GEN_PDV.hpp"

namespace moris::gen
{
    class PDV_Value : public PDV
    {

      private:
        real mValue;

      public:
        /**
         * Constructor
         *
         * @param aValue Constant value for this PDV
         */
        explicit PDV_Value( real aValue );

        /**
         * Gets if this PDV is active
         *
         * @return If PDV is active
         */
        bool is_active() override;

        /**
         * Get the PDV value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Coordinate values
         * @return Current value of this PDV
         */
        real get_value( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Get the PDV sensitivity with respect to ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Coordinate values
         * @return Vector of sensitivities to be returned
         */
        Matrix< DDRMat > get_sensitivities( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets the IDs of ADVs which this PDV depends on.
         *
         * @return ADV IDs
         */
        Vector< sint > get_determining_adv_ids( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates ) override;

    };
}

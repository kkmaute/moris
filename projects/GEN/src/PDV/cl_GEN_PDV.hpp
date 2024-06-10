/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

namespace moris::gen
{
    class PDV
    {
      private:
        moris_id mID       = gNoID;

      public:

        /**
         * constructor
         */
        PDV() = default;

      public:
        /**
         * trivial destructor
         */
        virtual ~PDV() = default;

        /**
         * Sets the PDV ID of this PDV
         *
         * @param aID
         */
        void set_id( moris_id aID );

        /**
         * Gets the PDV ID of this PDV
         *
         * @return
         */
        moris_id get_id();

        /**
         * Gets if this PDV is active
         *
         * @return If PDV is active
         */
        virtual bool is_active() = 0;

        /**
         * Get the PDV value
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Coordinate values
         * @return Current value of this PDV
         */
        virtual real get_value( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates ) = 0;

        /**
         * Get the PDV sensitivity with respect to ADVs
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Coordinate values
         * @return Vector of sensitivities to be returned
         */
        virtual Matrix< DDRMat > get_sensitivities( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates ) = 0;

        /**
         * Gets the IDs of ADVs which this PDV depends on.
         *
         * @return ADV IDs
         */
        virtual Vector< sint > get_determining_adv_ids( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates ) = 0;
    };
}

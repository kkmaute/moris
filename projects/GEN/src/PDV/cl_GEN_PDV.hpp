/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV.hpp
 *
 */

#ifndef MORIS_CL_GEN_PDV_HPP_
#define MORIS_CL_GEN_PDV_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class PDV
        {
          public:
            bool     mIsActive = true;
            moris_id mId       = gNoID;

          protected:
            /**
             * constructor
             */
            PDV();

          public:
            /**
             * trivial destructor
             */
            virtual ~PDV();

            /**
             * set PDV IDr
             */
            void
            set_id( const moris_id& tId )
            {
                mId = tId;
            };

            /**
             * get PDV ID
             */
            moris_id
            get_id()
            {
                return mId;
            };

            /**
             * Get the PDV value
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @return Current value of this PDV
             */
            virtual real get_value( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates ) = 0;

            virtual void
            set_value( const moris::real& aValue )
            {
                MORIS_ERROR( false, "PDV::set_value(), not implemented for this pdv type" );
            }

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
            virtual Matrix< DDSMat > get_determining_adv_ids( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates ) = 0;
        };
    }    // namespace ge
}    // namespace moris

#endif /* MORIS_CL_GEN_PDV_HPP_ */

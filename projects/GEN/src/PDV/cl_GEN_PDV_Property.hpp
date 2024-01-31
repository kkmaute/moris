/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV_Property.hpp
 *
 */

#ifndef MORIS_CL_GEN_PDV_PROPERTY_HPP
#define MORIS_CL_GEN_PDV_PROPERTY_HPP

#include "cl_GEN_PDV.hpp"
#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {
        class PDV_Property : public PDV
        {

          private:
            std::shared_ptr< Property > mProperty;

          public:
            /**
             * Constructor
             *
             * @param aPropertyPointer a GEN property pointer
             */
            PDV_Property( std::shared_ptr< Property > aPropertyPointer );

            /**
             * Destructor
             */
            ~PDV_Property();

            /**
             * Get the PDV value
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @return Current value of this PDV
             */
            real
            get_value(
                    uint                    aNodeIndex,
                    const Matrix< DDRMat >& aCoordinates );

            /**
             * Get the PDV sensitivity with respect to ADVs
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @return Vector of sensitivities to be returned
             */
            Matrix< DDRMat > get_sensitivities(
                    uint                    aNodeIndex,
                    const Matrix< DDRMat >& aCoordinates );

            /**
             * Gets the IDs of ADVs which this PDV depends on.
             *
             * @return ADV IDs
             */
            Matrix< DDSMat > get_determining_adv_ids(
                    uint                    aNodeIndex,
                    const Matrix< DDRMat >& aCoordinates );
        };
    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_GEN_PDV_PROPERTY_HPP

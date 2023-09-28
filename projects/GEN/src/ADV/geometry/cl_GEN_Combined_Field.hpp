/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Multifield.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"

namespace moris::ge
{
    class Combined_Field : public Field
    {

    private:
        Cell< std::shared_ptr< Field > > mFields;
        real mScale;

    public:

        /**
         * Combined field constructor
         *
         * @param aFields Created fields
         * @param aUseMinimum Whether or not to use minimum or maximum value when combining fields
         */
        Combined_Field(
                Cell< std::shared_ptr< Field > > aFields,
                bool                             aUseMinimum = true );

        /**
         * Given a node coordinate, the geometry needs to return the distance to the nearest function.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates vector of coordinate values
         * @return distance to nearest function
         */
        real get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node index or coordinate, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        const Matrix<DDRMat>& get_dfield_dadvs(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        virtual void get_dfield_dcoordinates(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities);

    };
}

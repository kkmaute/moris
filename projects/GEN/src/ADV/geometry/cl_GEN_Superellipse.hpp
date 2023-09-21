/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Superellipse.hpp
 *
 */

#pragma once

#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris::ge
{
    class Superellipse : public Field_Analytic
    {
    private:
        real mEpsilon = 1E-8;

    public:

        /**
         * Constructor, sets the pointers to advs and constant parameters for evaluations.
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADV vector
         * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aParameters Additional parameters
         */
        template <typename Vector_Type>
        Superellipse(
                Vector_Type&              aADVs,
                Matrix<DDUMat>            aGeometryVariableIndices,
                Matrix<DDUMat>            aADVIndices,
                Matrix<DDRMat>            aConstants,
              Level_Set_Parameters aParameters = {})
                : Field_Analytic( aADVs, aGeometryVariableIndices, aADVIndices, aConstants )
        {
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstants.length() == 8,
                        "A GEN Super-ellipse must be created with a total of exactly 8 variables (ADVs + constants).");

            MORIS_ERROR(*(mVariables(2)) > 0 and *(mVariables(3)) > 0,
                        "A GEN Super-ellipse must be created with positive semi-diameters.");

            MORIS_ERROR(std::abs(std::fmod(*(mVariables(4)),2.0)) < 1e-12,
                        "A GEN Super-ellipse must be created with an even exponent.");
        }

        /**
         * Constructor with only constant parameters
         *
         * @param aXCenter x-coordinate of the center of the superellipse
         * @param aYCenter y-coordiante of the center of the superellipse
         * @param aXSemidiameter Superellipse semi-diameter in the x direction
         * @param aYSemidiameter Superellipse semi-diameter in the y direction
         * @param aExponent Superellipse exponent
         * @param aParameters Additional parameters
         */
        Superellipse(
                real                      aXCenter,
                real                      aYCenter,
                real                      aXSemidiameter,
                real                      aYSemidiameter,
                real                      aExponent,
                real                      aScaling,
                real                      aRegularization,
                real                      aShift,
                Level_Set_Parameters aParameters = {});

        /**
         * Given a node coordinate, returns the field value.
         *
         * @param aCoordinates Coordinate values
         * @return Distance to this geometry
         */
        real get_field_value(const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node coordinate, evaluates the sensitivity of the geometry field with respect to all of the
         * geometry variables.
         *
         * @param aCoordinates Coordinate values
         * @return Vector of sensitivities
         */
        const Matrix<DDRMat>& get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates);

        /**
         * Given nodal coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities);

    };
}

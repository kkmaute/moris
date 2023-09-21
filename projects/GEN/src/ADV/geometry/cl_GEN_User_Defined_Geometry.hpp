/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_User_Defined_Geometry.hpp
 *
 */

#pragma once

#include "cl_GEN_User_Defined_Field.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"

namespace moris::ge
{
    class User_Defined_Geometry : public User_Defined_Field
    {
    public:

        /**
         * Constructor, sets the pointers to advs and constant parameters for evaluations.
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADV vector
         * @param aGeometryVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aFieldFunction User-defined function for evaluating the field field
         * @param aSensitivityFunction User-defined function for evaluating the field sensitivities
         * @param aParameters Additional parameters
         */
        template <typename Vector_Type>
        User_Defined_Geometry(
                Vector_Type&              aADVs,
                Matrix<DDUMat>            aGeometryVariableIndices,
                Matrix<DDUMat>            aADVIndices,
                Matrix<DDRMat>            aConstants,
                Field_Function            aFieldFunction,
                Sensitivity_Function      aSensitivityFunction,
                Level_Set_Parameters aParameters = {})
                : User_Defined_Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstants, aFieldFunction, aSensitivityFunction, aParameters)
        {
        }

        /**
         * Constructor with only constants and no sensitivities.
         *
         * @param aConstants The constant field variables not filled by ADVs
         * @param aFieldFunction User-defined function for evaluating the field field
         * @param aParameters Additional parameters
         */
        User_Defined_Geometry(
                Matrix<DDRMat>            aConstants,
                Field_Function            aFieldFunction,
                Level_Set_Parameters aParameters = {});

    };
}

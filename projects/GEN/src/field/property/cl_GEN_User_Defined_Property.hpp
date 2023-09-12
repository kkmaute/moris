/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_User_Defined_Property.hpp
 *
 */

#ifndef MORIS_CL_GEN_USER_DEFINED_PROPERTY_HPP
#define MORIS_CL_GEN_USER_DEFINED_PROPERTY_HPP

#include "cl_GEN_User_Defined_Field.hpp"
#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {
        class User_Defined_Property : public User_Defined_Field, public Property
        {
        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations.
             *
             * @tparam Vector_Type Type of vector where ADVs are stored
             * @param aADVs ADV vector
             * @param aPropertyVariableIndices Indices of field variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the field variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aFieldFunction User-defined function for evaluating the field field
             * @param aSensitivityFunction User-defined function for evaluating the field sensitivities
             * @param aParameters Additional parameters
             */
            template <typename Vector_Type>
            User_Defined_Property(
                    Vector_Type&              aADVs,
                    Matrix<DDUMat>            aPropertyVariableIndices,
                    Matrix<DDUMat>            aADVIndices,
                    Matrix<DDRMat>            aConstants,
                    Field_Function            aFieldFunction,
                    Sensitivity_Function      aSensitivityFunction,
                    Property_Field_Parameters aParameters = {})
                    : Field(aADVs, aPropertyVariableIndices, aADVIndices, aConstants, aParameters)
                    , User_Defined_Field(aADVs, aPropertyVariableIndices, aADVIndices, aConstants, aFieldFunction, aSensitivityFunction, aParameters)
                    , Property(aParameters)
            {
            }

            /**
             * Constructor with only constants and no sensitivities.
             *
             * @param aConstants The constant field variables not filled by ADVs
             * @param aFieldFunction User-defined function for evaluating the field field
             * @param aParameters Additional parameters
             */
            User_Defined_Property(
                    Matrix<DDRMat>            aConstants,
                    Field_Function            aFieldFunction,
                    Property_Field_Parameters aParameters = {});

        };
    }
}

#endif //MORIS_CL_GEN_USER_DEFINED_PROPERTY_HPP


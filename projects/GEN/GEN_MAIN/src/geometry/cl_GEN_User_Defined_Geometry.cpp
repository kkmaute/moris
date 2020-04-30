//
// Created by christopherson on 4/29/20.
//

#include "cl_GEN_User_Defined_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Geometry::User_Defined_Geometry(Matrix<DDRMat>& aADVs,
                                                     Matrix<DDUMat> aGeometryVariableIndices,
                                                     Matrix<DDUMat> aADVIndices,
                                                     Matrix<DDRMat> aConstantParameters,
                                                     MORIS_GEOMETRY_FUNCTION aFieldEvaluationFunction,
                                                     MORIS_GEOMETRY_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction)
                                                     : Geometry_Analytic(aADVs,
                                                             aGeometryVariableIndices,
                                                             aADVIndices,
                                                             aConstantParameters)
        {
            evaluate_field_value_user_defined = aFieldEvaluationFunction;
            evaluate_sensitivity_user_defined = aSensitivityEvaluationFunction;
        }

        //--------------------------------------------------------------------------------------------------------------

        real User_Defined_Geometry::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return this->evaluate_field_value_user_defined(aCoordinates, mGeometryVariables);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> User_Defined_Geometry::evaluate_sensitivity(const Matrix<DDRMat>& aCoordinates)
        {
            return this->evaluate_sensitivity_user_defined(aCoordinates, mGeometryVariables);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}



//
// Created by christopherson on 4/29/20.
//

#ifndef MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP
#define MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP

#include "cl_GEN_Geometry_Analytic.hpp"
#include "cl_Matrix.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        class User_Defined_Geometry : public Geometry_Analytic
        {

        private:
            MORIS_GEOMETRY_FUNCTION evaluate_field_value_user_defined;
            MORIS_GEOMETRY_SENSITIVITY_FUNCTION evaluate_sensitivity_user_defined;

        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             */
            User_Defined_Geometry(Matrix<DDRMat>& aADVs,
                                  Matrix<DDUMat> aGeometryVariableIndices,
                                  Matrix<DDUMat> aADVIndices,
                                  Matrix<DDRMat> aConstantParameters,
                                  MORIS_GEOMETRY_FUNCTION aFieldEvaluationFunction,
                                  MORIS_GEOMETRY_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction);

            /**
             * Given a node coordinate, the geometry needs to return the distance to the nearest function.
             *
             * @param aCoordinates vector of coordinate values
             * @return distance to nearest function
             */
            real evaluate_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate @param[in] aCoordinates, the function returns a matrix of relevant node coordinates
             * Where each row represents a design variable and each column is x, y, z sensitivities
             *
             * @param aCoordinates vector of coordinate values
             * @return matrix of sensitivities
             */
            Matrix<DDRMat> evaluate_sensitivity(const Matrix<DDRMat>& aCoordinates);

        };
    }
}

#endif //MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP

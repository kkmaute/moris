#ifndef MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP
#define MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        class User_Defined_Geometry : public Geometry, public Field_Analytic
        {

        private:
            MORIS_GEN_FIELD_FUNCTION evaluate_field_value_user_defined;
            MORIS_GEN_SENSITIVITY_FUNCTION evaluate_sensitivity_user_defined;

        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param aSensitivityEvaluationFunction User-defined function for evaluating the field sensitivities
             */
            User_Defined_Geometry(Matrix<DDRMat>& aADVs,
                                  Matrix<DDUMat> aGeometryVariableIndices,
                                  Matrix<DDUMat> aADVIndices,
                                  Matrix<DDRMat> aConstantParameters,
                                  MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                                  MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction = nullptr);

            /**
             * Constructor with only constant parameters
             *
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param aSensitivityEvaluationFunction User-defined function for evaluating the field sensitivities
             */
            User_Defined_Geometry(Matrix<DDRMat> aConstantParameters,
                                  MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                                  MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction = nullptr);

            /**
             * Given a node coordinate, the geometry needs to return the distance to the nearest function.
             *
             * @param aCoordinates vector of coordinate values
             * @return distance to nearest function
             */
            real evaluate_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate @param aCoordinates, the function returns a matrix of sensitivities of the
             * geometry location with respect to the ADVs
             *
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivities Matrix of sensitivities
             */
            void evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);

        private:

            /**
             * Sets the user-defined functions. Eliminates redundant code since it's the same logic for all constructors.
             *
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param aSensitivityEvaluationFunction User-defined function for evaluating the field sensitivities
             */
            void set_user_defined_functions(MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                                            MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction);

            /**
             * Used internally to automatically error out if no sensitivities were provided
             */
            static void no_sensitivities(const Matrix<DDRMat>&  aCoordinates,
                                         const Cell<real*>&     aParameters,
                                         Matrix<DDRMat>&        aSensitivities);

        };
    }
}

#endif //MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP

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
            MORIS_GEN_FIELD_FUNCTION get_field_value_user_defined;
            MORIS_GEN_SENSITIVITY_FUNCTION get_field_sensitivities_user_defined;

        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations.
             *
             * @tparam Vector_Type Type of vector where ADVs are stored
             * @param aADVs ADV vector
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param tSensitivitiesEvaluationFunction User-defined function for evaluating the field sensitivities
             * @param aParameters Additional parameters
             */
            template <typename Vector_Type>
            User_Defined_Geometry(
                    Vector_Type&                   aADVs,
                    Matrix<DDUMat>                 aGeometryVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstants,
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                    Field_Parameters               aParameters = {})
                    : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstants, aParameters)
            {
                this->set_user_defined_functions(aFieldEvaluationFunction, aSensitivityEvaluationFunction);
            }

            /**
             * Constructor with only constant parameters
             *
             * @param aConstants The constant field variables not filled by ADVs
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param aParameters Additional parameters
             */
            User_Defined_Geometry(
                    Matrix<DDRMat>           aConstants,
                    MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                    Field_Parameters         aParameters = {});

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
            const Matrix<DDRMat>& get_field_sensitivities(const Matrix<DDRMat>& aCoordinates);

        private:

            /**
             * Sets the user-defined functions. Eliminates redundant code since it's the same logic for all constructors.
             *
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param tSensitivitiesEvaluationFunction User-defined function for evaluating the field sensitivities
             */
            void set_user_defined_functions(
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction);

            /**
             * Used internally to automatically error out if no sensitivities were provided
             */
            static void no_sensitivities(
                    const Matrix<DDRMat>& aCoordinates,
                    const Cell<real*>&    aParameters,
                    Matrix<DDRMat>&       aSensitivities);

        };
    }
}

#endif //MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP

#ifndef MORIS_CL_GEN_SUPERELLIPSOID_HPP
#define MORIS_CL_GEN_SUPERELLIPSOID_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {
        class Superellipsoid : public Geometry, public Field_Analytic
        {
        private:
            real mEpsilon = 1E-8;

        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @tparam Vector_Type Type of vector where ADVs are stored
             * @param aADVs ADV vector
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            template <typename Vector_Type>
            Superellipsoid(
                    Vector_Type&              aADVs,
                    Matrix<DDUMat>            aGeometryVariableIndices,
                    Matrix<DDUMat>            aADVIndices,
                    Matrix<DDRMat>            aConstants,
                    Geometry_Field_Parameters aParameters = {})
                    : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstants, aParameters)
                    , Geometry(aParameters)
            {
                MORIS_ERROR(aGeometryVariableIndices.length() + aConstants.length() == 7,
                            "A GEN Superellipsoid must be created with a total of exactly 7 variables (ADVs + constants).");
                MORIS_ERROR(*(mFieldVariables(3)) > 0 and *(mFieldVariables(4)) > 0 and *(mFieldVariables(5)) > 0,
                            "A GEN Superellipsoid must be created with positive semidiameters.");
            }

            /**
             * Constructor with only constant parameters
             *
             * @param aXCenter x-coordinate of the center of the superellipsoid
             * @param aYCenter y-coordiante of the center of the superellipsoid
             * @param aZCenter z-coordinate of the center of the superellipsoid
             * @param aXSemidiameter Superellipsoid semi-diameter in the x direction
             * @param aYSemidiameter Superellipsoid semi-diameter in the y direction
             * @param aZSemidiameter Superellipsoid semi-diameter in the z direction
             * @param aExponent Superellipsoid exponent
             * @param aParameters Additional parameters
             */
            Superellipsoid(
                    real                      aXCenter,
                    real                      aYCenter,
                    real                      aZCenter,
                    real                      aXSemidiameter,
                    real                      aYSemidiameter,
                    real                      aZSemidiameter,
                    real                      aExponent,
                    Geometry_Field_Parameters aParameters = {});

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

        };
    }
}

#endif //MORIS_CL_GEN_SUPERELLIPSOID_HPP

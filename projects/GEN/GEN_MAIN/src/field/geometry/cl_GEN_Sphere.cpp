#include "cl_GEN_Sphere.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Sphere::Sphere(
                Matrix<DDRMat>&  aADVs,
                Matrix<DDUMat>   aGeometryVariableIndices,
                Matrix<DDUMat>   aADVIndices,
                Matrix<DDRMat>   aConstants,
                Field_Parameters aParameters)
                : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstants, aParameters)
        {
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstants.length() == 4,
                        "A GEN Sphere must be created with a total of exactly 4 variables (ADVs + constant parameters)");
        }

        //--------------------------------------------------------------------------------------------------------------

        Sphere::Sphere(
                sol::Dist_Vector* aOwnedADVs,
                Matrix<DDUMat>    aGeometryVariableIndices,
                Matrix<DDUMat>    aADVIndices,
                Matrix<DDRMat>    aConstants,
                Field_Parameters  aParameters)
                : Field(aOwnedADVs, aGeometryVariableIndices, aADVIndices, aConstants, aParameters)
        {
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstants.length() == 4,
                        "A GEN Sphere must be created with a total of exactly 4 variables (ADVs + constant parameters)");
        }

        //--------------------------------------------------------------------------------------------------------------

        Sphere::Sphere(real aXCenter, real aYCenter, real aZCenter, real aRadius, Field_Parameters aParameters)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aRadius}}), aParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Sphere::get_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tZCenter = *(mFieldVariables(2));
            real tRadius = *(mFieldVariables(3));

            // Evaluate field
            return sqrt(pow(aCoordinates(0) - tXCenter, 2)
                      + pow(aCoordinates(1) - tYCenter, 2)
                      + pow(aCoordinates(2) - tZCenter, 2)) - tRadius;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Sphere::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tZCenter = *(mFieldVariables(2));

            // Calculate sensitivities
            real tConstant = sqrt(pow(aCoordinates(0) - tXCenter, 2)
                    + pow(aCoordinates(1) - tYCenter, 2)
                    + pow(aCoordinates(2) - tZCenter, 2));
            tConstant = tConstant ? 1 / tConstant : 0.0;
            mSensitivities(0) = tConstant * (tXCenter - aCoordinates(0));
            mSensitivities(1) = tConstant * (tYCenter - aCoordinates(1));
            mSensitivities(2) = tConstant * (tZCenter - aCoordinates(2));
            mSensitivities(3) = -1;

            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

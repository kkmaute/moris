#include "cl_GEN_Superellipsoid.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Superellipsoid::Superellipsoid(
                Matrix<DDRMat>& aADVs,
                Matrix<DDUMat>  aGeometryVariableIndices,
                Matrix<DDUMat>  aADVIndices,
                Matrix<DDRMat>  aConstantParameters,
                std::string     aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aNumPatterns,
                sint            aRefinementFunctionIndex,
                sint            aBSplineMeshIndex,
                real            aBSplineLowerBound,
                real            aBSplineUpperBound)
                : Field(aADVs,
                        aGeometryVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 7,
                        "A GEN Superellipsoid must be created with a total of exactly 7 variables (ADVs + constant parameters).");
            MORIS_ERROR(*(mFieldVariables(3)) > 0 and *(mFieldVariables(4)) > 0 and *(mFieldVariables(5)) > 0,
                        "A GEN Superellipsoid must be created with positive semidiameters.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Superellipsoid::Superellipsoid(
                sol::Dist_Vector* aOwnedADVs,
                Matrix<DDUMat>    aGeometryVariableIndices,
                Matrix<DDUMat>    aADVIndices,
                Matrix<DDRMat>    aConstantParameters,
                std::string       aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aNumPatterns,
                sint              aRefinementFunctionIndex,
                sint              aBSplineMeshIndex,
                real              aBSplineLowerBound,
                real              aBSplineUpperBound)
                : Field(aOwnedADVs,
                        aGeometryVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 7,
                        "A GEN Superellipsoid must be created with a total of exactly 7 variables (ADVs + constant parameters).");
            MORIS_ERROR(*(mFieldVariables(3)) > 0 and *(mFieldVariables(4)) > 0 and *(mFieldVariables(5)) > 0,
                        "A GEN Superellipsoid must be created with positive semidiameters.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Superellipsoid::Superellipsoid(
                real        aXCenter,
                real        aYCenter,
                real        aZCenter,
                real        aXSemidiameter,
                real        aYSemidiameter,
                real        aZSemidiameter,
                real        aExponent,
                std::string aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aNumPatterns,
                sint        aRefinementFunctionIndex,
                sint        aBSplineMeshIndex,
                real        aBSplineLowerBound,
                real        aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aXSemidiameter, aYSemidiameter, aZSemidiameter, aExponent}}),
                        aName,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            MORIS_ERROR(*(mFieldVariables(3)) > 0 and *(mFieldVariables(4)) > 0 and *(mFieldVariables(5)) > 0,
                        "A GEN Superellipsoid must be created with positive semidiameters.");
        }

        //--------------------------------------------------------------------------------------------------------------

        real Superellipsoid::get_field_value_geometry(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tZCenter = *(mFieldVariables(2));
            real tXSemidiameter = *(mFieldVariables(3));
            real tYSemidiameter = *(mFieldVariables(4));
            real tZSemidiameter = *(mFieldVariables(5));
            real tExponent = *(mFieldVariables(6));

            // Evaluate field
            return pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent)
                     + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent)
                     + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent), 1.0 / tExponent) - 1.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Superellipsoid::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tZCenter = *(mFieldVariables(2));
            real tXSemidiameter = *(mFieldVariables(3));
            real tYSemidiameter = *(mFieldVariables(4));
            real tZSemidiameter = *(mFieldVariables(5));
            real tExponent = *(mFieldVariables(6));

            // Constant in all calculations
            real tConstant = pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent)
                           + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent)
                           + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent);
            tConstant = tConstant ? pow(tConstant, -1.0 + (1.0 / tExponent)) : 0.0;

            // Calculate sensitivities
            Matrix<DDRMat> tSensitivities(1, 7);
            tSensitivities(0) = -tConstant * pow(1.0 / tXSemidiameter, tExponent) * (aCoordinates(0) - tXCenter)
                    * pow(std::abs(tXCenter - aCoordinates(0)), tExponent - 2.0);
            tSensitivities(1) = -tConstant * pow(1.0 / tYSemidiameter, tExponent) * (aCoordinates(1) - tYCenter)
                    * pow(std::abs(tYCenter - aCoordinates(1)), tExponent - 2.0);
            tSensitivities(2) = -tConstant * pow(1.0 / tZSemidiameter, tExponent) * (aCoordinates(2) - tZCenter)
                    * pow(std::abs(tZCenter - aCoordinates(2)), tExponent - 2.0);
            tSensitivities(3) = -tConstant * pow(1.0 / tXSemidiameter, tExponent + 1.0)
                    * pow(std::abs(tXCenter - aCoordinates(0)), tExponent);
            tSensitivities(4) = -tConstant * pow(1.0 / tYSemidiameter, tExponent + 1.0)
                    * pow(std::abs(tYCenter - aCoordinates(1)), tExponent);
            tSensitivities(5) = -tConstant * pow(1.0 / tZSemidiameter, tExponent + 1.0)
                    * pow(std::abs(tZCenter - aCoordinates(2)), tExponent);

            // TODO? this uses FD only because the analytical function is super complicated for the exponent derivative
            tSensitivities(6) = (pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent + mEpsilon)
                                   + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent + mEpsilon)
                                   + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent + mEpsilon), 1.0 / (tExponent + mEpsilon))
                               - pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent - mEpsilon)
                                   + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent - mEpsilon)
                                   + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent - mEpsilon), 1.0 / (tExponent - mEpsilon)))
                               / (2.0 * mEpsilon);

            return tSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


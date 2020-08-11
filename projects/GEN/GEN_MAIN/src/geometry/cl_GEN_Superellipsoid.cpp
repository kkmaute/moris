#include "cl_GEN_Superellipsoid.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Superellipsoid::Superellipsoid(Matrix<DDRMat>& aADVs,
                                       Matrix<DDUMat>  aGeometryVariableIndices,
                                       Matrix<DDUMat>  aADVIndices,
                                       Matrix<DDRMat>  aConstantParameters,
                                       sint            aNumRefinements,
                                       sint            aRefinementFunctionIndex,
                                       sint            aBSplineMeshIndex,
                                       real            aBSplineLowerBound,
                                       real            aBSplineUpperBound)
                : Field(aADVs,
                        aGeometryVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aNumRefinements,
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

        Superellipsoid::Superellipsoid(real aXCenter,
                                       real aYCenter,
                                       real aZCenter,
                                       real aXSemidiameter,
                                       real aYSemidiameter,
                                       real aZSemidiameter,
                                       real aExponent,
                                       sint aNumRefinements,
                                       sint aRefinementFunctionIndex,
                                       sint aBSplineMeshIndex,
                                       real aBSplineLowerBound,
                                       real aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aXSemidiameter, aYSemidiameter, aZSemidiameter, aExponent}}),
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            MORIS_ERROR(*(mFieldVariables(3)) > 0 and *(mFieldVariables(4)) > 0 and *(mFieldVariables(5)) > 0,
                        "A GEN Superellipsoid must be created with positive semidiameters.");
        }

        //--------------------------------------------------------------------------------------------------------------

        real Superellipsoid::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
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

        void Superellipsoid::evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
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
            aSensitivities.set_size(1, 7);
            aSensitivities(0) = -tConstant * pow(1.0 / tXSemidiameter, tExponent) * (aCoordinates(0) - tXCenter)
                    * pow(std::abs(tXCenter - aCoordinates(0)), tExponent - 2.0);
            aSensitivities(1) = -tConstant * pow(1.0 / tYSemidiameter, tExponent) * (aCoordinates(1) - tYCenter)
                    * pow(std::abs(tYCenter - aCoordinates(1)), tExponent - 2.0);
            aSensitivities(2) = -tConstant * pow(1.0 / tZSemidiameter, tExponent) * (aCoordinates(2) - tZCenter)
                    * pow(std::abs(tZCenter - aCoordinates(2)), tExponent - 2.0);
            aSensitivities(3) = -tConstant * pow(1.0 / tXSemidiameter, tExponent + 1.0)
                    * pow(std::abs(tXCenter - aCoordinates(0)), tExponent);
            aSensitivities(4) = -tConstant * pow(1.0 / tYSemidiameter, tExponent + 1.0)
                    * pow(std::abs(tYCenter - aCoordinates(1)), tExponent);
            aSensitivities(5) = -tConstant * pow(1.0 / tZSemidiameter, tExponent + 1.0)
                    * pow(std::abs(tZCenter - aCoordinates(2)), tExponent);

            // TODO? this uses FD only because the analytical function is super complicated for the exponent derivative
            aSensitivities(6) = (pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent + mEpsilon)
                                   + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent + mEpsilon)
                                   + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent + mEpsilon), 1.0 / (tExponent + mEpsilon))
                               - pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent - mEpsilon)
                                   + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent - mEpsilon)
                                   + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent - mEpsilon), 1.0 / (tExponent - mEpsilon)))
                               / (2.0 * mEpsilon);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


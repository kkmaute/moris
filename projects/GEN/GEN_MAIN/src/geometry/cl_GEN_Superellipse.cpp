#include "cl_GEN_Superellipse.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Superellipse::Superellipse(Matrix<DDRMat>& aADVs,
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
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 5,
                        "A GEN Superellipse must be created with a total of exactly 5 variables (ADVs + constant parameters).");
            MORIS_ERROR(*(mFieldVariables(2)) > 0 and *(mFieldVariables(3)) > 0,
                        "A GEN Superellipse must be created with positive semidiameters.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Superellipse::Superellipse(real aXCenter,
                                   real aYCenter,
                                   real aXSemidiameter,
                                   real aYSemidiameter,
                                   real aExponent,
                                   sint aNumRefinements,
                                   sint aRefinementFunctionIndex,
                                   sint aBSplineMeshIndex,
                                   real aBSplineLowerBound,
                                   real aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aXSemidiameter, aYSemidiameter, aExponent}}),
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            MORIS_ERROR(*(mFieldVariables(2)) > 0 and *(mFieldVariables(3)) > 0,
                        "A GEN Superellipse must be created with positive semidiameters.");
        }

        //--------------------------------------------------------------------------------------------------------------

        real Superellipse::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tXSemidiameter = *(mFieldVariables(2));
            real tYSemidiameter = *(mFieldVariables(3));
            real tExponent = *(mFieldVariables(4));

            // Evaluate field
            return pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent)
                     + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent), 1.0 / tExponent) - 1.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Superellipse::evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tXSemidiameter = *(mFieldVariables(2));
            real tYSemidiameter = *(mFieldVariables(3));
            real tExponent = *(mFieldVariables(4));

            // Constant in all calculations
            real tConstant = pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent)
                           + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent);
            tConstant = tConstant ? pow(tConstant, -1.0 + (1.0 / tExponent)) : 0.0;

            // Calculate sensitivities
            aSensitivities.set_size(1, 5);
            aSensitivities(0) = -tConstant * pow(1.0 / tXSemidiameter, tExponent) * (aCoordinates(0) - tXCenter)
                    * pow(std::abs(tXCenter - aCoordinates(0)), tExponent - 2.0);
            aSensitivities(1) = -tConstant * pow(1.0 / tYSemidiameter, tExponent) * (aCoordinates(1) - tYCenter)
                    * pow(std::abs(tYCenter - aCoordinates(1)), tExponent - 2.0);
            aSensitivities(2) = -tConstant * pow(1.0 / tXSemidiameter, tExponent + 1.0)
                    * pow(std::abs(tXCenter - aCoordinates(0)), tExponent);
            aSensitivities(3) = -tConstant * pow(1.0 / tYSemidiameter, tExponent + 1.0)
                    * pow(std::abs(tYCenter - aCoordinates(1)), tExponent);

            // TODO? this uses FD only because the analytical function is super complicated for the exponent derivative
            aSensitivities(4) = (pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent + mEpsilon)
                                   + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent + mEpsilon), 1.0 / (tExponent + mEpsilon))
                               - pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent - mEpsilon)
                                   + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent - mEpsilon), 1.0 / (tExponent - mEpsilon)))
                               / (2.0 * mEpsilon);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

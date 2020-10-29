#include "cl_GEN_Circle.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Circle::Circle(Matrix<DDRMat>& aADVs,
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
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 3,
                        "A GEN Circle must be created with a total of exactly 3 variables (ADVs + constant parameters)");
        }

        //--------------------------------------------------------------------------------------------------------------

        Circle::Circle(sol::Dist_Vector* aOwnedADVs,
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
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 3,
                        "A GEN Circle must be created with a total of exactly 3 variables (ADVs + constant parameters)");
        }

        //--------------------------------------------------------------------------------------------------------------

        Circle::Circle(real        aXCenter,
                       real        aYCenter,
                       real        aRadius,
                       std::string       aName,
                       Matrix<DDSMat>  aNumRefinements,
                       Matrix<DDSMat>  aNumPatterns,
                       sint        aRefinementFunctionIndex,
                       sint        aBSplineMeshIndex,
                       real        aBSplineLowerBound,
                       real        aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aRadius}}),
                        aName,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Circle::get_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tRadius = *(mFieldVariables(2));

            // Evaluate field
            return sqrt(pow(aCoordinates(0) - tXCenter, 2) + pow(aCoordinates(1) - tYCenter, 2)) - tRadius;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Circle::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            // Initialize
            Matrix<DDRMat> tSensitivities(1, 3);

            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));

            // Calculate sensitivities
            real tConstant = sqrt(pow(aCoordinates(0) - tXCenter, 2) + pow(aCoordinates(1) - tYCenter, 2));
            tConstant = tConstant ? 1 / tConstant : 0.0;
            tSensitivities(0) = tConstant * (tXCenter - aCoordinates(0));
            tSensitivities(1) = tConstant * (tYCenter - aCoordinates(1));
            tSensitivities(2) = -1.0;

            return tSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

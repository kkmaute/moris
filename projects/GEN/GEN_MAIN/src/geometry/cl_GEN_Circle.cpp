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
                       sint            aNumRefinements,
                       sint            aRefinementFunctionIndex,
                       sint            aBSplineMeshIndex,
                       real            aLevelSetLowerBound,
                       real            aLevelSetUpperBound)
                : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstantParameters),
                  Geometry(aNumRefinements,
                           aRefinementFunctionIndex,
                           aBSplineMeshIndex,
                           aLevelSetLowerBound,
                           aLevelSetUpperBound)
        {
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 3,
                        "A circle geometry must be created with a total of exactly 3 adv and constant parameters");
        }

        //--------------------------------------------------------------------------------------------------------------

        Circle::Circle(real aXCenter,
                       real aYCenter,
                       real aRadius,
                       sint aNumRefinements,
                       sint aRefinementFunctionIndex,
                       sint aBSplineMeshIndex,
                       real aLevelSetLowerBound,
                       real aLevelSetUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aRadius}})),
                  Geometry(aNumRefinements,
                           aRefinementFunctionIndex,
                           aBSplineMeshIndex,
                           aLevelSetLowerBound,
                           aLevelSetUpperBound)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Circle::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tRadius = *(mFieldVariables(2));

            Matrix<DDRMat> tCenter(2, 1);
            tCenter(0) = tXCenter;
            tCenter(1) = tYCenter;

            // Evaluate field
            return sqrt(pow(aCoordinates(0) - tXCenter, 2) + pow(aCoordinates(1) - tYCenter, 2)) - tRadius;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Circle::evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));

            // Calculate sensitivities
            aSensitivities.resize(1, 3);
            real tDenominator = sqrt(pow(aCoordinates(0) - tXCenter, 2) + pow(aCoordinates(1) - tYCenter, 2));
            if (tDenominator == 0.0)
            {
                aSensitivities(0) = 0.0;
                aSensitivities(1) = 0.0;
            }
            else
            {
                aSensitivities(0) = (tXCenter - aCoordinates(0)) / tDenominator;
                aSensitivities(1) = (tYCenter - aCoordinates(1)) / tDenominator;
            }
            aSensitivities(2) = -1.0;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

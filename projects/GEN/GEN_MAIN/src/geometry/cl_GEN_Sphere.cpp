#include "cl_GEN_Sphere.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Sphere::Sphere(Matrix<DDRMat>& aADVs,
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
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 4,
                        "A sphere geometry must be created with a total of exactly 4 adv and constant parameters");
        }

        //--------------------------------------------------------------------------------------------------------------

        Sphere::Sphere(real aXCenter,
                       real aYCenter,
                       real aZCenter,
                       real aRadius,
                       sint aNumRefinements,
                       sint aRefinementFunctionIndex,
                       sint aBSplineMeshIndex,
                       real aLevelSetLowerBound,
                       real aLevelSetUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aRadius}})),
                  Geometry(aNumRefinements,
                           aRefinementFunctionIndex,
                           aBSplineMeshIndex,
                           aLevelSetLowerBound,
                           aLevelSetUpperBound)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Sphere::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
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

        void Sphere::evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tZCenter = *(mFieldVariables(2));

            // Calculate sensitivities
            aSensitivities.resize(1, 4);
            real tDenominator = sqrt(pow(aCoordinates(0) - tXCenter, 2)
                                   + pow(aCoordinates(1) - tYCenter, 2)
                                   + pow(aCoordinates(2) - tZCenter, 2));
            if (tDenominator == 0.0)
            {
                aSensitivities(0) = 0.0;
                aSensitivities(1) = 0.0;
                aSensitivities(2) = 0.0;
            }
            else
            {
                aSensitivities(0) = (tXCenter - aCoordinates(0)) / tDenominator;
                aSensitivities(1) = (tYCenter - aCoordinates(1)) / tDenominator;
                aSensitivities(2) = (tZCenter - aCoordinates(2)) / tDenominator;
            }
            aSensitivities(3) = -1;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

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
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 4,
                        "A GEN Sphere must be created with a total of exactly 4 variables (ADVs + constant parameters)");
        }

        //--------------------------------------------------------------------------------------------------------------

        Sphere::Sphere(sol::Dist_Vector* aOwnedADVs,
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
            MORIS_ERROR(aGeometryVariableIndices.length() + aConstantParameters.length() == 4,
                        "A GEN Sphere must be created with a total of exactly 4 variables (ADVs + constant parameters)");
        }

        //--------------------------------------------------------------------------------------------------------------

        Sphere::Sphere(real        aXCenter,
                       real        aYCenter,
                       real        aZCenter,
                       real        aRadius,
                       std::string aName,
                       Matrix<DDSMat>  aNumRefinements,
                       Matrix<DDSMat>  aNumPatterns,
                       sint        aRefinementFunctionIndex,
                       sint        aBSplineMeshIndex,
                       real        aBSplineLowerBound,
                       real        aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aRadius}}),
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

        Matrix<DDRMat> Sphere::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tZCenter = *(mFieldVariables(2));

            // Calculate sensitivities
            Matrix<DDRMat> tSensitivities(1, 4);
            real tConstant = sqrt(pow(aCoordinates(0) - tXCenter, 2)
                    + pow(aCoordinates(1) - tYCenter, 2)
                    + pow(aCoordinates(2) - tZCenter, 2));
            tConstant = tConstant ? 1 / tConstant : 0.0;
            tSensitivities(0) = tConstant * (tXCenter - aCoordinates(0));
            tSensitivities(1) = tConstant * (tYCenter - aCoordinates(1));
            tSensitivities(2) = tConstant * (tZCenter - aCoordinates(2));
            tSensitivities(3) = -1;

            return tSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

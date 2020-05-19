#include "cl_GEN_Circle.hpp"

#include "fn_norm.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Circle::Circle(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters)
        : Geometry_Analytic(aADVs, aGeometryVariableIndices, aADVIndices, aConstantParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Circle::Circle(real aXCenter, real aYCenter, real aRadius) : Geometry_Analytic(Matrix<DDRMat>({{aXCenter, aYCenter, aRadius}}))
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Circle::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            Matrix<DDRMat> tCenter(1, 2);
            tCenter(0) = *(mFieldVariables(0));
            tCenter(1) = *(mFieldVariables(1));

            // Evaluate field
            moris::real tFunctionValue = norm(aCoordinates - tCenter) - *(mFieldVariables(2));
            
            return tFunctionValue;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Circle::evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            // Initialize sensitivity matrix
            aSensitivities.resize(3, 2);

            // Get variables
            moris::real tXCenter = *(mFieldVariables(0));
            moris::real tYCenter = *(mFieldVariables(1));
            moris::real tRadius = *(mFieldVariables(2));

            // dx/dr
            // Set sign based on value under square root
            moris::real sign = 0.0;
            moris::real tSqrt = tRadius * tRadius - std::pow((aCoordinates(1) - tYCenter), 2);
            if(tSqrt < 0.0)
            {
                sign = -1.0;
            }
            else if(tSqrt > 0.0)
            {
                sign = 1.0;
            }

            // Calculate
            aSensitivities(0, 0) = sign * tRadius / std::sqrt(std::abs(tSqrt));

            // dy/dr
            // Set sign based on value under square root
            tSqrt = tRadius * tRadius - std::pow((aCoordinates(0) - tXCenter),2);
            if(tSqrt < 0.0)
            {
                sign = -1.0;
            }
            else if(tSqrt > 0.0)
            {
                sign = 1.0;
            }

            // Calculate
            aSensitivities(0, 1) = sign * tRadius / std::sqrt(std::abs(tSqrt));

            // Fill remaining values in tSensitivity
            aSensitivities(1,0) = 1.0; // dx/dxc
            aSensitivities(1,1) = 0.0; // dy/dxc
            aSensitivities(2,0) = 0.0; // dx/dyc
            aSensitivities(2,1) = 1.0; // dy/dyc
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


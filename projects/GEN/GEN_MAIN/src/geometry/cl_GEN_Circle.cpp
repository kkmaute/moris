#include "cl_GEN_Circle.hpp"

#include "fn_norm.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        real Circle::evaluate_field_value(const moris::Matrix<moris::DDRMat>& aCoordinates)
        {
            // Get variables
            Matrix<DDRMat> tCenter(2, 1);
            tCenter(0) = *(mGeometryVariables(0));
            tCenter(1) = *(mGeometryVariables(1));

            // Evaluate field
            moris::real tFunctionValue = norm( aCoordinates - tCenter ) - *(mGeometryVariables(2));
            
            return tFunctionValue;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat> Circle::evaluate_sensitivity(const moris::Matrix<moris::DDRMat>& aCoordinates)
        {
            // Initialize sensitivity matrix
            moris::Matrix< moris::DDRMat > tSensitivityDxDp(3, 2, 0.0);

            // Get variables
            moris::real tXCenter = *(mGeometryVariables(0));
            moris::real tYCenter = *(mGeometryVariables(1));
            moris::real tRadius = *(mGeometryVariables(2));

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
            tSensitivityDxDp(0, 0) = sign * tRadius / std::sqrt(std::abs(tSqrt));

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
            tSensitivityDxDp(0, 1) = sign * tRadius / std::sqrt(std::abs(tSqrt));

            // Fill remaining values in tSensitivity
            tSensitivityDxDp(1,0) = 1.0; // dx/dxc
            tSensitivityDxDp(1,1) = 0.0; // dy/dxc
            tSensitivityDxDp(2,0) = 0.0; // dx/dyc
            tSensitivityDxDp(2,1) = 1.0; // dy/dyc

            return tSensitivityDxDp;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


//
// Created by christopherson on 4/17/20.
//

#include "cl_GEN_Plane.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------
        
        Plane::Plane(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters)
        : Geometry_Analytic(aADVs, aGeometryVariableIndices, aADVIndices, aConstantParameters)
        {
            if (mGeometryVariables.size() == 4)
            {
                m_eval_field = &Plane::eval_field_2d;
                m_eval_sensitivity = &Plane::eval_sensitivity_2d;
            }
            else if (mGeometryVariables.size() == 6)
            {
                m_eval_field = &Plane::eval_field_3d;
                m_eval_sensitivity = &Plane::eval_sensitivity_3d;
            }
            else
            {
                MORIS_ERROR(false, "Incorrect number of parameters passed for construction of a ge::Plane");
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Plane::Plane(real aXCenter, real aYCenter, real aZCenter, real aXNormal, real aYNormal, real aZNormal)
        : Geometry_Analytic(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aXNormal, aYNormal, aZNormal}})) 
        {
            m_eval_field = &Plane::eval_field_3d;
            m_eval_sensitivity = &Plane::eval_sensitivity_3d;
        }

        //--------------------------------------------------------------------------------------------------------------
    
        Plane::Plane(real aXCenter, real aYCenter, real aXNormal, real aYNormal)
        : Geometry_Analytic(Matrix<DDRMat>({{aXCenter, aYCenter, aXNormal, aYNormal}}))
        {
            m_eval_field = &Plane::eval_field_2d;
            m_eval_sensitivity = &Plane::eval_sensitivity_2d;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Plane::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return (this->*m_eval_field)(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Plane::evaluate_sensitivity(const Matrix<DDRMat>& aCoordinates)
        {
            return (this->*m_eval_sensitivity)(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::real Plane::eval_field_2d(Matrix<DDRMat> const & aCoordinates)
        {
            real tXCenter = *(mGeometryVariables(0));
            real tYCenter = *(mGeometryVariables(1));
            real tXNormal = *(mGeometryVariables(2));
            real tYNormal = *(mGeometryVariables(3));

            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter);
        }

        //--------------------------------------------------------------------------------------------------------------
        
        moris::real Plane::eval_field_3d(Matrix<DDRMat> const & aCoordinates)
        {
            real tXCenter = *(mGeometryVariables(0));
            real tYCenter = *(mGeometryVariables(1));
            real tZCenter = *(mGeometryVariables(2));
            real tXNormal = *(mGeometryVariables(3));
            real tYNormal = *(mGeometryVariables(4));
            real tZNormal = *(mGeometryVariables(5));

            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter) + tZNormal * (aCoordinates(2) - tZCenter);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Plane::eval_sensitivity_2d(Matrix<DDRMat> const & aCoordinates)
        {
            MORIS_ERROR(false, "Sensitivities not implemented for 2d plane.");
            return Matrix<DDRMat>(1, 1, 0.0);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Plane::eval_sensitivity_3d(Matrix<DDRMat> const & aCoordinates)
        {
            MORIS_ERROR(false, "Sensitivities not implemented for 3d plane.");
            return Matrix<DDRMat>(1, 1, 0.0);
        }

        //--------------------------------------------------------------------------------------------------------------
        
    }
}

//
// Created by christopherson on 4/17/20.
//

#include "cl_GEN_Plane.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Plane::Plane(Matrix<DDRMat>& aADVs,
                     Matrix<DDUMat>  aGeometryVariableIndices,
                     Matrix<DDUMat>  aADVIndices,
                     Matrix<DDRMat>  aConstantParameters,
                     sint            aNumRefinements,
                     sint            aRefinementFunctionIndex)
                : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstantParameters),
                  Geometry(aNumRefinements, aRefinementFunctionIndex)
        {
            if (mFieldVariables.size() == 4)
            {
                m_eval_field = &Plane::eval_field_2d;
                m_eval_sensitivity = &Plane::eval_sensitivity_2d;
            }
            else if (mFieldVariables.size() == 6)
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

        Plane::Plane(real aXCenter,
                     real aYCenter,
                     real aZCenter,
                     real aXNormal,
                     real aYNormal,
                     real aZNormal,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aXNormal, aYNormal, aZNormal}})),
                  Geometry(aNumRefinements, aRefinementFunctionIndex)
        {
            m_eval_field = &Plane::eval_field_3d;
            m_eval_sensitivity = &Plane::eval_sensitivity_3d;
        }

        //--------------------------------------------------------------------------------------------------------------
    
        Plane::Plane(real aXCenter,
                     real aYCenter,
                     real aXNormal,
                     real aYNormal,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aXNormal, aYNormal}})),
                  Geometry(aNumRefinements, aRefinementFunctionIndex)
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

        void Plane::evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            (this->*m_eval_sensitivity)(aCoordinates, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::real Plane::eval_field_2d(Matrix<DDRMat> const & aCoordinates)
        {
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tXNormal = *(mFieldVariables(2));
            real tYNormal = *(mFieldVariables(3));

            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter);
        }

        //--------------------------------------------------------------------------------------------------------------
        
        moris::real Plane::eval_field_3d(Matrix<DDRMat> const & aCoordinates)
        {
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tZCenter = *(mFieldVariables(2));
            real tXNormal = *(mFieldVariables(3));
            real tYNormal = *(mFieldVariables(4));
            real tZNormal = *(mFieldVariables(5));

            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter) + tZNormal * (aCoordinates(2) - tZCenter);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Plane::eval_sensitivity_2d(Matrix<DDRMat> const & aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            MORIS_ERROR(false, "Sensitivities not implemented for 2d plane.");
        }

        //--------------------------------------------------------------------------------------------------------------

        void Plane::eval_sensitivity_3d(Matrix<DDRMat> const & aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            MORIS_ERROR(false, "Sensitivities not implemented for 3d plane.");
        }

        //--------------------------------------------------------------------------------------------------------------
        
    }
}

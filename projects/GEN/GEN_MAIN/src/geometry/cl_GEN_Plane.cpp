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
                MORIS_ERROR(false, "Incorrect number of parameters passed for construction of a GEN Plane.");
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Plane::Plane(sol::Dist_Vector* aOwnedADVs,
                     Matrix<DDUMat>    aGeometryVariableIndices,
                     Matrix<DDUMat>    aADVIndices,
                     Matrix<DDRMat>    aConstantParameters,
                     sint              aNumRefinements,
                     sint              aRefinementFunctionIndex,
                     sint              aBSplineMeshIndex,
                     real              aBSplineLowerBound,
                     real              aBSplineUpperBound)
                : Field(aOwnedADVs,
                        aGeometryVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
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
                MORIS_ERROR(false, "Incorrect number of parameters passed for construction of a GEN Plane.");
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
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aXNormal, aYNormal, aZNormal}}),
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
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
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aXNormal, aYNormal}}),
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            m_eval_field = &Plane::eval_field_2d;
            m_eval_sensitivity = &Plane::eval_sensitivity_2d;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Plane::get_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return (this->*m_eval_field)(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Plane::evaluate_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            (this->*m_eval_sensitivity)(aCoordinates, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::real Plane::eval_field_2d(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tXNormal = *(mFieldVariables(2));
            real tYNormal = *(mFieldVariables(3));

            // Evaluate field value
            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter);
        }

        //--------------------------------------------------------------------------------------------------------------
        
        moris::real Plane::eval_field_3d(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tZCenter = *(mFieldVariables(2));
            real tXNormal = *(mFieldVariables(3));
            real tYNormal = *(mFieldVariables(4));
            real tZNormal = *(mFieldVariables(5));

            // Evaluate field value
            return tXNormal * (aCoordinates(0) - tXCenter) + tYNormal * (aCoordinates(1) - tYCenter) + tZNormal * (aCoordinates(2) - tZCenter);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Plane::eval_sensitivity_2d(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tXNormal = *(mFieldVariables(2));
            real tYNormal = *(mFieldVariables(3));

            // Evaluate sensitivities
            aSensitivities.set_size(1, 4);
            aSensitivities(0) = -tXNormal;
            aSensitivities(1) = -tYNormal;
            aSensitivities(2) = aCoordinates(0) - tXCenter;
            aSensitivities(3) = aCoordinates(1) - tYCenter;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Plane::eval_sensitivity_3d(Matrix<DDRMat> const & aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            MORIS_ERROR(false, "Sensitivities not implemented for 3d plane.");
        }

        //--------------------------------------------------------------------------------------------------------------
        
    }
}

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
                     std::string       aName,
                     Matrix<DDSMat>    aNumRefinements,
                     Matrix<DDSMat>    aNumPatterns,
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

        Plane::Plane(real           aXCenter,
                     real           aYCenter,
                     real           aZCenter,
                     real           aXNormal,
                     real           aYNormal,
                     real           aZNormal,
                     std::string    aName,
                     Matrix<DDSMat> aNumRefinements,
                     Matrix<DDSMat> aNumPatterns,
                     sint           aRefinementFunctionIndex,
                     sint           aBSplineMeshIndex,
                     real           aBSplineLowerBound,
                     real           aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aXNormal, aYNormal, aZNormal}}),
                        aName,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            m_eval_field = &Plane::eval_field_3d;
            m_eval_sensitivity = &Plane::eval_sensitivity_3d;
        }

        //--------------------------------------------------------------------------------------------------------------
    
        Plane::Plane(real           aXCenter,
                     real           aYCenter,
                     real           aXNormal,
                     real           aYNormal,
                     std::string    aName,
                     Matrix<DDSMat> aNumRefinements,
                     Matrix<DDSMat> aNumPatterns,
                     sint           aRefinementFunctionIndex,
                     sint           aBSplineMeshIndex,
                     real           aBSplineLowerBound,
                     real           aBSplineUpperBound)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aXNormal, aYNormal}}),
                        aName,
                        aNumRefinements,
                        aNumPatterns,
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

        Matrix<DDRMat> Plane::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            return (this->*m_eval_sensitivity)(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        real Plane::eval_field_2d(const Matrix<DDRMat>& aCoordinates)
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
        
        real Plane::eval_field_3d(const Matrix<DDRMat>& aCoordinates)
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

        Matrix<DDRMat> Plane::eval_sensitivity_2d(const Matrix<DDRMat>& aCoordinates)
        {
            // Get variables
            real tXCenter = *(mFieldVariables(0));
            real tYCenter = *(mFieldVariables(1));
            real tXNormal = *(mFieldVariables(2));
            real tYNormal = *(mFieldVariables(3));

            // Evaluate sensitivities
            Matrix<DDRMat> tSensitivities(1, 4);
            tSensitivities(0) = -tXNormal;
            tSensitivities(1) = -tYNormal;
            tSensitivities(2) = aCoordinates(0) - tXCenter;
            tSensitivities(3) = aCoordinates(1) - tYCenter;

            return tSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Plane::eval_sensitivity_3d(const Matrix<DDRMat>& aCoordinates)
        {
            MORIS_ERROR(false, "Sensitivities not implemented for 3d plane.");
            return {{}};
        }

        //--------------------------------------------------------------------------------------------------------------
        
    }
}

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(Matrix<DDRMat>& aADVs,
                     Matrix<DDUMat> aFieldVariableIndices,
                     Matrix<DDUMat> aADVIndices,
                     Matrix<DDRMat> aConstantParameters,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : mFieldVariables(aFieldVariableIndices.length() + aConstantParameters.length()),
                  mConstantParameters(aConstantParameters),
                  mADVDependencies(aFieldVariableIndices.length() + aConstantParameters.length(), 1, -1),
                  mDependsOnADVs(aADVIndices.length()),
                  mNumADVs(aADVs.length()),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)

        {
            // Resize field variables
            this->assign_adv_dependencies(aFieldVariableIndices, aADVIndices);

            // Fill with pointers to ADVs
            for (uint tADVFillIndex = 0; tADVFillIndex < aFieldVariableIndices.length(); tADVFillIndex++)
            {
                mFieldVariables(aFieldVariableIndices(tADVFillIndex)) = &(aADVs(aADVIndices(tADVFillIndex)));
            }
            
            // Fill constant parameters
            this->fill_constant_parameters();
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(sol::Dist_Vector* aOwnedADVs,
                     Matrix<DDUMat> aFieldVariableIndices,
                     Matrix<DDUMat> aADVIndices,
                     Matrix<DDRMat> aConstantParameters,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : mFieldVariables(aFieldVariableIndices.length() + aConstantParameters.length()),
                  mConstantParameters(aConstantParameters),
                  mADVDependencies(aFieldVariableIndices.length() + aConstantParameters.length(), 1, -1),
                  mDependsOnADVs(aADVIndices.length()),
                  mNumADVs(aOwnedADVs->vec_local_length()),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)

        {
            // Resize field variables
            this->assign_adv_dependencies(aFieldVariableIndices, aADVIndices);

            // Fill with pointers to ADVs
            for (uint tADVFillIndex = 0; tADVFillIndex < aFieldVariableIndices.length(); tADVFillIndex++)
            {
                mFieldVariables(aFieldVariableIndices(tADVFillIndex)) = &(*aOwnedADVs)(aADVIndices(tADVFillIndex));
            }

            // Fill constant parameters
            this->fill_constant_parameters();
        }

        //--------------------------------------------------------------------------------------------------------------
        
        Field::Field(Matrix<DDRMat>& aADVs,
                     uint aStartingADVIndex,
                     uint aNumFieldVariables,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : mFieldVariables(aNumFieldVariables),
                  mADVDependencies(aNumFieldVariables, 1),
                  mDependsOnADVs(true),
                  mNumADVs(aADVs.length()),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)
        {
            // Check for ADV size
            MORIS_ERROR((aStartingADVIndex + aNumFieldVariables) <= aADVs.length(),
                        "GEN field constructor with number of field variables given can only be called with an "
                        "ADV vector that has already been resized to an adequate length.");

            // Set variables from ADVs
            for (uint tVariableIndex = 0; tVariableIndex < aNumFieldVariables; tVariableIndex++)
            {
                mFieldVariables(tVariableIndex) = &(aADVs(aStartingADVIndex + tVariableIndex));
                mADVDependencies(tVariableIndex) = aStartingADVIndex + tVariableIndex;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(sol::Dist_Vector* aOwnedADVs,
                     uint aStartingADVIndex,
                     uint aNumFieldVariables,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : mFieldVariables(aNumFieldVariables),
                  mADVDependencies(aNumFieldVariables, 1),
                  mDependsOnADVs(true),
                  mNumADVs(aOwnedADVs->vec_local_length()),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)
        {
            // Check for ADV size
            MORIS_ERROR((aStartingADVIndex + aNumFieldVariables) <= (uint)aOwnedADVs->vec_local_length(),
                        "GEN field constructor with number of field variables given can only be called with an "
                        "ADV distributed vector with adequate local length.");

            // Set variables from ADVs
            for (uint tVariableIndex = 0; tVariableIndex < aNumFieldVariables; tVariableIndex++)
            {
                mFieldVariables(tVariableIndex) = &(*aOwnedADVs)(aStartingADVIndex + tVariableIndex + 1);
                mADVDependencies(tVariableIndex) = aStartingADVIndex + tVariableIndex;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(Matrix<DDRMat> aConstantParameters,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : mFieldVariables(aConstantParameters.length()),
                  mConstantParameters(aConstantParameters),
                  mADVDependencies(aConstantParameters.length(), 1, -1),
                  mDependsOnADVs(false),
                  mNumADVs(0),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)
        {
            // Fill field variables with all constant parameters
            for (uint tVariableIndex = 0; tVariableIndex < mConstantParameters.length(); tVariableIndex++)
            {
                mFieldVariables(tVariableIndex) = &(mConstantParameters(tVariableIndex));
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field()
        {
            MORIS_ERROR(false, "The default constructor of a GEN Field should never be used. It only exists because "
                               "the compiler on Blanca is stupid.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::~Field()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::evaluate_sensitivity(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities)
        {
            // Resize sensitivities
            aSensitivities.set_size(1, mNumADVs, 0.0);

            // Evaluate all sensitivities
            Matrix<DDRMat> tTempSensitivities;
            this->evaluate_all_sensitivities(aNodeIndex, aCoordinates, tTempSensitivities);

            // Return only what is needed
            for (uint tSensitivityIndex = 0; tSensitivityIndex < mADVDependencies.length(); tSensitivityIndex++)
            {
                if (mADVDependencies(tSensitivityIndex) >= 0)
                {
                    aSensitivities(mADVDependencies(tSensitivityIndex)) = tTempSensitivities(tSensitivityIndex);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::reset_child_nodes()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Field::depends_on_advs()
        {
            return mDependsOnADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Field::get_num_refinements()
        {
            return mNumRefinements;
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Field::get_refinement_function_index()
        {
            return mRefinementFunctionIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Field::get_bspline_mesh_index()
        {
            return mBSplineMeshIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field::get_bspline_lower_bound()
        {
            return mBSplineLowerBound;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field::get_bspline_upper_bound()
        {
            return mBSplineUpperBound;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Field::conversion_to_bsplines()
        {
            return (mBSplineMeshIndex >= 0);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::assign_adv_dependencies(
                Matrix<DDUMat>& aFieldVariableIndices,
                Matrix<DDUMat>& aADVIndices)
        {
            // Check that the number of field variables indices equals the number of ADV indices
            MORIS_ERROR(aFieldVariableIndices.length() == aADVIndices.length(),
                        "Number of field variables indices must equal the number of ADV indices in a GEN field.");

            // Set ADV dependencies
            for (uint tADVFillIndex = 0; tADVFillIndex < aFieldVariableIndices.length(); tADVFillIndex++)
            {
                mADVDependencies(aFieldVariableIndices(tADVFillIndex)) = aADVIndices(tADVFillIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Field::fill_constant_parameters()
        {
            uint tParameterIndex = 0;
            for (uint tVariableIndex = 0; tVariableIndex < mFieldVariables.size(); tVariableIndex++)
            {
                if (mFieldVariables(tVariableIndex) == nullptr)
                {
                    mFieldVariables(tVariableIndex) = &(mConstantParameters(tParameterIndex++));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


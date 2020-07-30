#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(Matrix<DDRMat>& aADVs,
                     Matrix<DDUMat> aFieldVariableIndices,
                     Matrix<DDUMat> aADVIndices,
                     Matrix<DDRMat> aConstantParameters)
                : mADVIndices(aADVIndices),
                  mConstantParameters(aConstantParameters),
                  mActiveVariables(aFieldVariableIndices.length() + mConstantParameters.length(), true)

        {
            // Check that the number of field variables indices equals the number of ADV indices, resize field variables
            MORIS_ERROR(aFieldVariableIndices.length() == aADVIndices.length(),
                        "Number of field variables indices must equal the number of ADV indices in a GEN field.");

            // Store number of ADVs
            mNumADVs = aADVs.length();

            // Resize field variables
            uint tNumInputs = aFieldVariableIndices.length() + mConstantParameters.length();
            mFieldVariables.resize(tNumInputs);

            // Fill with pointers to ADVs
            for (uint tADVFillIndex = 0; tADVFillIndex < aFieldVariableIndices.length(); tADVFillIndex++)
            {
                mFieldVariables(aFieldVariableIndices(tADVFillIndex)) = &(aADVs(aADVIndices(tADVFillIndex)));
            }

            // Fill with constant parameters and identify these variables
            uint tParameterIndex = 0;
            for (uint tVariableIndex = 0; tVariableIndex < tNumInputs; tVariableIndex++)
            {
                if (mFieldVariables(tVariableIndex) == nullptr)
                {
                    mFieldVariables(tVariableIndex) = &(mConstantParameters(tParameterIndex++));
                    mActiveVariables(tVariableIndex) = false;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        Field::Field(Matrix<DDRMat>& aADVs,
                     uint aADVIndex,
                     uint aNumFieldVariables)
                : mActiveVariables(aNumFieldVariables, true),
                  mNumADVs(aADVs.length())
        {
            // Check for ADV size
            MORIS_ERROR((aADVIndex + aNumFieldVariables) <= aADVs.length(),
                        "GEN field constructor with number of field variables given can only be called with an "
                        "ADV vector that has already been resized to an adequate length.");

            // Resize ADVs/variables
            mADVIndices.resize(aNumFieldVariables, 1);
            mFieldVariables.resize(aNumFieldVariables);

            // Set variables from ADVs
            for (uint tVariableIndex = 0; tVariableIndex < aNumFieldVariables; tVariableIndex++)
            {
                mADVIndices(tVariableIndex) = aADVIndex + tVariableIndex;
                mFieldVariables(tVariableIndex) = &(aADVs(aADVIndex + tVariableIndex));
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(Matrix<DDRMat> aConstantParameters)
                : mConstantParameters(aConstantParameters),
                  mActiveVariables(aConstantParameters.length(), false),
                  mNumADVs(0)
        {
            // Resize field variables
            uint tNumInputs = mConstantParameters.length();
            mFieldVariables.resize(tNumInputs);

            // Fill field variables
            for (uint tVariableIndex = 0; tVariableIndex < tNumInputs; tVariableIndex++)
            {
                mFieldVariables(tVariableIndex) = &(mConstantParameters(tVariableIndex));
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::~Field()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::evaluate_sensitivity(      uint            aIndex,
                                         const Matrix<DDRMat>& aCoordinates,
                                               Matrix<DDRMat>& aSensitivities)
        {
            // Resize sensitivities
            aSensitivities.set_size(1, mNumADVs, 0.0);

            // Evaluate all sensitivities
            Matrix<DDRMat> tTempSensitivities;
            this->evaluate_all_sensitivities(aIndex, aCoordinates, tTempSensitivities);

            // Return only what is needed
            for (uint tSensitivityIndex = 0; tSensitivityIndex < mActiveVariables.size(); tSensitivityIndex++)
            {
                if (mActiveVariables(tSensitivityIndex))
                {
                    aSensitivities(mADVIndices(tSensitivityIndex)) = tTempSensitivities(tSensitivityIndex);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Field::depends_on_advs()
        {
            return mADVIndices.length();
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


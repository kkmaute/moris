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
                : mConstantParameters(aConstantParameters),
                  mFieldADVDependencies(aFieldVariableIndices.length() + mConstantParameters.length(), 1, -1),
                  mDependsOnADVs(aADVIndices.length()),
                  mNumADVs(aADVs.length()),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)

        {
            // Check that the number of field variables indices equals the number of ADV indices, resize field variables
            MORIS_ERROR(aFieldVariableIndices.length() == aADVIndices.length(),
                        "Number of field variables indices must equal the number of ADV indices in a GEN field.");

            // Resize field variables
            uint tNumInputs = aFieldVariableIndices.length() + mConstantParameters.length();
            mFieldVariables.resize(tNumInputs);

            // Fill with pointers to ADVs
            for (uint tADVFillIndex = 0; tADVFillIndex < aFieldVariableIndices.length(); tADVFillIndex++)
            {
                mFieldVariables(aFieldVariableIndices(tADVFillIndex)) = &(aADVs(aADVIndices(tADVFillIndex)));
                mFieldADVDependencies(aFieldVariableIndices(tADVFillIndex)) = aADVIndices(tADVFillIndex);
            }

            // Fill with constant parameters and identify these variables
            uint tParameterIndex = 0;
            for (uint tVariableIndex = 0; tVariableIndex < tNumInputs; tVariableIndex++)
            {
                if (mFieldVariables(tVariableIndex) == nullptr)
                {
                    mFieldVariables(tVariableIndex) = &(mConstantParameters(tParameterIndex++));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        Field::Field(Matrix<DDRMat>& aADVs,
                     uint aADVIndex,
                     uint aNumFieldVariables,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : mFieldVariables(aNumFieldVariables),
                  mFieldADVDependencies(aNumFieldVariables, 1),
                  mDependsOnADVs(true),
                  mNumADVs(aADVs.length()),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)
        {
            // Check for ADV size
            MORIS_ERROR((aADVIndex + aNumFieldVariables) <= aADVs.length(),
                        "GEN field constructor with number of field variables given can only be called with an "
                        "ADV vector that has already been resized to an adequate length.");

            // Set variables from ADVs
            for (uint tVariableIndex = 0; tVariableIndex < aNumFieldVariables; tVariableIndex++)
            {
                mFieldADVDependencies(tVariableIndex) = aADVIndex + tVariableIndex;
                mFieldVariables(tVariableIndex) = &(aADVs(aADVIndex + tVariableIndex));
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(Matrix<DDRMat> aConstantParameters,
                     sint aNumRefinements,
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aBSplineLowerBound,
                     real aBSplineUpperBound)
                : mConstantParameters(aConstantParameters),
                  mFieldADVDependencies(aConstantParameters.length(), 1, -1),
                  mDependsOnADVs(false),
                  mNumADVs(0),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)
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

        void Field::evaluate_sensitivity(uint                  aNodeIndex,
                                         const Matrix<DDRMat>& aCoordinates,
                                         Matrix<DDRMat>&       aSensitivities)
        {
            // Resize sensitivities
            aSensitivities.set_size(1, mNumADVs, 0.0);

            // Evaluate all sensitivities
            Matrix<DDRMat> tTempSensitivities;
            this->evaluate_all_sensitivities(aNodeIndex, aCoordinates, tTempSensitivities);

            // Return only what is needed
            for (uint tSensitivityIndex = 0; tSensitivityIndex < mFieldADVDependencies.length(); tSensitivityIndex++)
            {
                if (mFieldADVDependencies(tSensitivityIndex) >= 0)
                {
                    aSensitivities(mFieldADVDependencies(tSensitivityIndex)) = tTempSensitivities(tSensitivityIndex);
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

    }
}


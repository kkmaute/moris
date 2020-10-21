#include "cl_GEN_Field.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(Matrix<DDRMat>& aADVs,
                     Matrix<DDUMat>  aFieldVariableIndices,
                     Matrix<DDUMat>  aADVIndices,
                     Matrix<DDRMat>  aConstantParameters,
                     sint            aNumRefinements,
                     sint            aRefinementFunctionIndex,
                     sint            aBSplineMeshIndex,
                     real            aBSplineLowerBound,
                     real            aBSplineUpperBound)
                : mFieldVariables(aFieldVariableIndices.length() + aConstantParameters.length()),
                  mConstantParameters(aConstantParameters),
                  mDeterminingADVIds(aFieldVariableIndices.length() + aConstantParameters.length(), 1, -1),
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
                     Matrix<DDUMat>    aFieldVariableIndices,
                     Matrix<DDUMat>    aADVIndices,
                     Matrix<DDRMat>    aConstantParameters,
                     sint              aNumRefinements,
                     sint              aRefinementFunctionIndex,
                     sint              aBSplineMeshIndex,
                     real              aBSplineLowerBound,
                     real              aBSplineUpperBound)
                : mFieldVariables(aFieldVariableIndices.length() + aConstantParameters.length()),
                  mConstantParameters(aConstantParameters),
                  mDeterminingADVIds(aFieldVariableIndices.length() + aConstantParameters.length(), 1, -1),
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

        Field::Field(const Matrix<DDSMat>& aSharedADVIds,
                     sint                  aNumRefinements,
                     sint                  aRefinementFunctionIndex,
                     sint                  aBSplineMeshIndex,
                     real                  aBSplineLowerBound,
                     real                  aBSplineUpperBound)
                : mFieldVariables(aSharedADVIds.length()),
                  mDeterminingADVIds(aSharedADVIds),
                  mDependsOnADVs(true),
                  mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mBSplineLowerBound(aBSplineLowerBound),
                  mBSplineUpperBound(aBSplineUpperBound)
        {
            // Create shared distributed vector
            sol::Matrix_Vector_Factory tDistributedFactory;
            std::shared_ptr<sol::Dist_Map> tSharedADVMap = tDistributedFactory.create_map(aSharedADVIds);
            mSharedADVs = tDistributedFactory.create_vector(tSharedADVMap);

            // Set variables from ADVs
            uint tNumSharedADVs = aSharedADVIds.length();
            for (uint tVariableIndex = 0; tVariableIndex < tNumSharedADVs; tVariableIndex++)
            {
                mFieldVariables(tVariableIndex) = &(*mSharedADVs)(aSharedADVIds(tVariableIndex));
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(Matrix<DDRMat> aConstantParameters,
                     sint           aNumRefinements,
                     sint           aRefinementFunctionIndex,
                     sint           aBSplineMeshIndex,
                     real           aBSplineLowerBound,
                     real           aBSplineUpperBound)
                : mFieldVariables(aConstantParameters.length()),
                  mConstantParameters(aConstantParameters),
                  mDeterminingADVIds(aConstantParameters.length(), 1, -1),
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

        void Field::import_advs(sol::Dist_Vector* aOwnedADVs)
        {
            if (mSharedADVs)
            {
                mSharedADVs->import_local_to_global(*aOwnedADVs);
                mNumADVs = aOwnedADVs->vec_local_length();
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

        bool Field::conversion_to_bsplines()
        {
            return (mBSplineMeshIndex >= 0);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Field::get_determining_adv_ids(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return mDeterminingADVIds;
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
                mDeterminingADVIds(aFieldVariableIndices(tADVFillIndex)) = aADVIndices(tADVFillIndex);
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


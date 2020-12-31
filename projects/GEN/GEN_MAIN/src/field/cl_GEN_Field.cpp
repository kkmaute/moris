#include "cl_GEN_Field.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------
        // Get addresses inside of different vector types
        //--------------------------------------------------------------------------------------------------------------

        real* get_address(Matrix<DDRMat>& aVector, uint aIndex)
        {
            return &aVector(aIndex);
        }

        real* get_address(sol::Dist_Vector* aVector, uint aIndex)
        {
            return &(*aVector)(aIndex);
        }

        //--------------------------------------------------------------------------------------------------------------
        // Definitions
        //--------------------------------------------------------------------------------------------------------------

        template <typename Vector_Type>
        Field::Field(Vector_Type&     aADVs,
                     Matrix<DDUMat>   aFieldVariableIndices,
                     Matrix<DDUMat>   aADVIndices,
                     Matrix<DDRMat>   aConstants,
                     Field_Parameters aParameters)
                : mFieldVariables(aFieldVariableIndices.length() + aConstants.length())
                , mSensitivities(1, aFieldVariableIndices.length() + aConstants.length())
                , mConstants(aConstants)
                , mParameters(aParameters)
                , mDeterminingADVIds(aFieldVariableIndices.length() + aConstants.length(), 1, -1)
                , mDependsOnADVs(aADVIndices.length())
        {
            // Check that refinement information is correct
            MORIS_ERROR(mParameters.mNumRefinements.length() == mParameters.mRefinementMeshIndices.length(),
                    "The entries given for number of refinements must line up with the number of refinement patterns.");

            // Check and assign ADV dependencies
            this->assign_adv_dependencies(aFieldVariableIndices, aADVIndices);

            // Fill with pointers to ADVs
            for (uint tADVFillIndex = 0; tADVFillIndex < aFieldVariableIndices.length(); tADVFillIndex++)
            {
                mFieldVariables(aFieldVariableIndices(tADVFillIndex)) = get_address(aADVs, aADVIndices(tADVFillIndex));
            }
            
            // Fill constant parameters
            this->fill_constant_parameters();
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(const Matrix<DDSMat>&  aSharedADVIds,
                     std::shared_ptr<Field> aField)
                : mFieldVariables(aSharedADVIds.length())
                , mSensitivities(1, aSharedADVIds.length())
                , mParameters(aField->mParameters)
                , mDeterminingADVIds(aSharedADVIds)
                , mDependsOnADVs(true)
        {
            // Check that refinement information is correct
            MORIS_ERROR(mParameters.mNumRefinements.length() == mParameters.mRefinementMeshIndices.length(),
                    "The entries given for number of refinements must line up with the number of refinement patterns.");

            // Create shared distributed vector
            sol::Matrix_Vector_Factory tDistributedFactory;
            sol::Dist_Map* tSharedADVMap = tDistributedFactory.create_map(aSharedADVIds);
            mSharedADVs = tDistributedFactory.create_vector(tSharedADVMap, 1, false, true);

            // Set variables from ADVs
            uint tNumSharedADVs = aSharedADVIds.length();
            for (uint tVariableIndex = 0; tVariableIndex < tNumSharedADVs; tVariableIndex++)
            {
                mFieldVariables(tVariableIndex) = &(*mSharedADVs)(aSharedADVIds(tVariableIndex));
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(Matrix<DDRMat>   aConstants,
                     Field_Parameters aParameters)
                : mFieldVariables(aConstants.length())
                , mSensitivities(1, aConstants.length())
                , mConstants(aConstants)
                , mParameters(aParameters)
                , mDeterminingADVIds(aConstants.length(), 1, -1)
                , mDependsOnADVs(false)
        {
            // Fill constant parameters
            this->fill_constant_parameters();
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(std::shared_ptr<Field> aField)
                : mParameters(aField->mParameters)
                , mDependsOnADVs(aField->mDependsOnADVs)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field()
        {
            MORIS_ERROR(false, "The default constructor of a GEN Field should never be used. It only exists because "
                               "the compiler on Blanca can't figure out that this constructor is not needed.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::~Field()
        {
            delete mSharedADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::import_advs(sol::Dist_Vector* aOwnedADVs)
        {
            if (mSharedADVs)
            {
                mSharedADVs->import_local_to_global(*aOwnedADVs);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::reset_nodal_information()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Field::store_field_values()
        {
            return (mParameters.mBSplineMeshIndex > -2);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Field::conversion_to_bsplines()
        {
            return (mParameters.mBSplineMeshIndex > -1);
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

        std::string Field::get_name()
        {
            return mParameters.mName;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDSMat > & Field::get_num_refinements()
        {
            return mParameters.mNumRefinements;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDSMat > & Field::get_refinement_mesh_indices()
        {
            return mParameters.mRefinementMeshIndices;
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Field::get_refinement_function_index()
        {
            return mParameters.mRefinementFunctionIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Field::get_bspline_mesh_index()
        {
            return (uint)mParameters.mBSplineMeshIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field::get_bspline_lower_bound()
        {
            return mParameters.mBSplineLowerBound;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field::get_bspline_upper_bound()
        {
            return mParameters.mBSplineUpperBound;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::assign_adv_dependencies(
                Matrix<DDUMat> aFieldVariableIndices,
                Matrix<DDUMat> aADVIndices)
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
                    mFieldVariables(tVariableIndex) = &(mConstants(tParameterIndex++));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        // Explicit template instantiation
        //--------------------------------------------------------------------------------------------------------------

        template
        Field::Field(Matrix<DDRMat>&  aADVs,
                     Matrix<DDUMat>   aFieldVariableIndices,
                     Matrix<DDUMat>   aADVIndices,
                     Matrix<DDRMat>   aConstants,
                     Field_Parameters aParameters);

        template
        Field::Field(sol::Dist_Vector*& aADVs,
                     Matrix<DDUMat>     aFieldVariableIndices,
                     Matrix<DDUMat>     aADVIndices,
                     Matrix<DDRMat>     aConstants,
                     Field_Parameters   aParameters);

        //--------------------------------------------------------------------------------------------------------------

    }
}


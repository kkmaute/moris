/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field.cpp
 *
 */

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
            this->set_advs(aADVs);

            // Fill constant parameters
            this->fill_constant_parameters();

            mLabel = mParameters.mName;
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(const Matrix<DDUMat>&  aFieldVariableIndices,
                     const Matrix<DDSMat>&  aSharedADVIds,
                     mtk::Mesh_Pair         aMeshPair,
                     std::shared_ptr<Field> aField)
                : mtk::Field(aMeshPair)
                , mFieldVariables(aSharedADVIds.length())
                , mSensitivities(1, aSharedADVIds.length())
                , mParameters(aField->mParameters)
                , mDeterminingADVIds(aSharedADVIds)
                , mDependsOnADVs(true)
        {
            // Check that refinement information is correct
            MORIS_ERROR(mParameters.mNumRefinements.length() == mParameters.mRefinementMeshIndices.length(),
                    "The entries given for number of refinements must line up with the number of refinement patterns.");

            // Check that the field variable indices match the shared ADV Ids
            MORIS_ERROR(aFieldVariableIndices.length() == aSharedADVIds.length(),
                    "Number of field variable indices must equal the number of ADV IDs in a GEN Field.");

            // Create shared distributed vector
            sol::Matrix_Vector_Factory tDistributedFactory;
            sol::Dist_Map* tSharedADVMap = tDistributedFactory.create_map(aSharedADVIds);
            mSharedADVs = tDistributedFactory.create_vector(tSharedADVMap, 1, false, true);

            // Set variables from ADVs
            uint tNumSharedADVs = aSharedADVIds.length();
            for (uint tVariable = 0; tVariable < tNumSharedADVs; tVariable++)
            {
                mFieldVariables(aFieldVariableIndices(tVariable)) = &(*mSharedADVs)(aSharedADVIds(tVariable));
            }

            mLabel = mParameters.mName;
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

            mLabel = mParameters.mName;
        }

        //--------------------------------------------------------------------------------------------------------------

        Field::Field(std::shared_ptr<Field> aField)
                : mParameters(aField->mParameters)
                , mDependsOnADVs(aField->mDependsOnADVs)
        {
            mLabel = mParameters.mName;
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

        template <typename Vector_Type>
        void Field::set_advs(Vector_Type& aADVs)
        {
            for (uint tVariableIndex = 0; tVariableIndex < mDeterminingADVIds.length(); tVariableIndex++)
            {
                if (mDeterminingADVIds(tVariableIndex) > -1)
                {
                    mFieldVariables(tVariableIndex) = get_address(aADVs, mDeterminingADVIds(tVariableIndex));
                }
            }
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

        void Field::add_nodal_data(mtk::Interpolation_Mesh* aMesh)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::reset_nodal_data()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Field::intended_storage()
        {
            return (mParameters.mDiscretizationMeshIndex > -2);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Field::intended_discretization()
        {
            return (mParameters.mDiscretizationMeshIndex > -1);
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

        moris_index Field::get_discretization_mesh_index() const
        {
            MORIS_ASSERT(mParameters.mDiscretizationMeshIndex >= 0,
                    "A discretization is not intended for this field. Check this with intended_discretization() first.");
            return mParameters.mDiscretizationMeshIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Field::get_discretization_order() const
        {
            return mMeshPair.get_interpolation_mesh()->get_discretization_order( this->get_discretization_mesh_index() );
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field::get_discretization_lower_bound()
        {
            return mParameters.mDiscretizationLowerBound;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field::get_discretization_upper_bound()
        {
            return mParameters.mDiscretizationUpperBound;
        }

        //--------------------------------------------------------------------------------------------------------------

        mtk::Mesh_Pair Field::get_mesh_pair()
        {
            return mtk::Mesh_Pair(nullptr, nullptr);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field::set_num_original_nodes( uint aNumOriginalNodes )
        {
            mNumOriginalNodes = aNumOriginalNodes;
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

        template
        void Field::set_advs(Matrix<DDRMat>& aADVs);

        template
        void Field::set_advs(sol::Dist_Vector*& aADVs);

        //--------------------------------------------------------------------------------------------------------------

    }
}


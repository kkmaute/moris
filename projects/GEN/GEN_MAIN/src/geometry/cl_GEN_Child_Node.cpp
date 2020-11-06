#include "cl_GEN_Child_Node.hpp"
#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Child_Node::Child_Node(Matrix<DDUMat>             aParentNodeIndices,
                               Cell<Matrix<DDRMat>>       aParentNodeCoordinates,
                               const xtk::Basis_Function& aBasisFunction,
                               Matrix<DDRMat>             aLocalCoordinates)
                : mParentNodeIndices(aParentNodeIndices),
                  mParentNodeCoordinates(aParentNodeCoordinates),
                  mLocalCoordinates(aLocalCoordinates)
        {
            // Shift local coordinates
            real tEpsilon = 1E-12;
            for (uint tDimension = 0; tDimension < mLocalCoordinates.length(); tDimension++)
            {
                mLocalCoordinates(tDimension) = std::min(mLocalCoordinates(tDimension), 1.0 - tEpsilon);
                mLocalCoordinates(tDimension) = std::max(mLocalCoordinates(tDimension), tEpsilon - 1.0);
            }

            // Evaluate basis function
            aBasisFunction.evaluate_basis_function(mLocalCoordinates, mBasisValues);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Child_Node::get_local_coordinates()
        {
            return mLocalCoordinates;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Child_Node::interpolate_field_value(Field* aField)
        {
            // Get field values from parents
            Matrix<DDRMat> tGeometryFieldValues(mParentNodeIndices.length(), 1);
            for (uint tParentNode = 0; tParentNode < mParentNodeIndices.length(); tParentNode++)
            {
                tGeometryFieldValues(tParentNode) = aField->get_field_value(mParentNodeIndices(tParentNode), mParentNodeCoordinates(tParentNode));
            }

            // Return interpolated value
            return Matrix<DDRMat>(mBasisValues * tGeometryFieldValues)(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Child_Node::join_field_sensitivities(Field* aField)
        {
            // Initialize using first parent
            mJoinedSensitivities = aField->get_field_sensitivities(
                    mParentNodeIndices(0),
                    mParentNodeCoordinates(0));
            mJoinedSensitivities = mJoinedSensitivities * mBasisValues(0);

            // Get sensitivity values from other parents
            for (uint tParentNode = 1; tParentNode < mParentNodeIndices.length(); tParentNode++)
            {
                // Get scaled sensitivities
                const Matrix<DDRMat>& tParentSensitivities = mBasisValues(tParentNode) * aField->get_field_sensitivities(
                        mParentNodeIndices(tParentNode),
                        mParentNodeCoordinates(tParentNode));
                
                // Join sensitivities
                uint tJoinedSensitivityLength = mJoinedSensitivities.n_cols();
                mJoinedSensitivities.resize(1, tJoinedSensitivityLength + tParentSensitivities.n_cols());
                for (uint tParentSensitivity = 0; tParentSensitivity < tParentSensitivities.n_cols(); tParentSensitivity++)
                {
                    mJoinedSensitivities(tJoinedSensitivityLength + tParentSensitivity) = tParentSensitivities(tParentSensitivity);
                }
            }

            return mJoinedSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Child_Node::join_determining_adv_ids(Field* aField)
        {
            // Initialize using first parent
            Matrix<DDSMat> tJoinedDeterminingADVs = aField->get_determining_adv_ids(
                    mParentNodeIndices(0),
                    mParentNodeCoordinates(0));

            // Get sensitivity values from other parents
            for (uint tParentNode = 1; tParentNode < mParentNodeIndices.length(); tParentNode++)
            {
                // Get scaled sensitivities
                Matrix<DDSMat> tParentDependingADVs = aField->get_determining_adv_ids(
                        mParentNodeIndices(tParentNode),
                        mParentNodeCoordinates(tParentNode));

                // Join sensitivities
                uint tJoinedADVLength = tJoinedDeterminingADVs.n_cols();
                tJoinedDeterminingADVs.resize(1, tJoinedADVLength + tParentDependingADVs.n_cols());
                for (uint tParentADV = 0; tParentADV < tParentDependingADVs.n_cols(); tParentADV++)
                {
                    tJoinedDeterminingADVs(tJoinedADVLength + tParentADV) = tParentDependingADVs(tParentADV);
                }
            }

            return tJoinedDeterminingADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
